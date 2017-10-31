#!/usr/bin/env python3
import os, sys, re
import struct
import gzip
from itertools import chain

# This is a stand-alone script to inspect a CBCL file - see the format description at
# https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2_guide_15051736_v2.pdf
# I'll get this working then use it as the basis for the tile reader update.

def main(cbcl_file):

    def p(key, val): print("  {}: {}".format(key, val))
    def w(msg): print("!! {} !!".format(msg))

    # First job is to open the file and inspect the header...
    with open(cbcl_file, 'rb') as fh:

        # First 12 bytes are fixed fields.
        header_bytes = fh.read(12)
        h_version, h_size, h_basebits, h_qbits, h_bins = struct.unpack('<HIBBI', header_bytes)

        print("Header info for {}...".format(cbcl_file))
        p("version", h_version)
        if h_version != 1: w("expected version 1")
        p("header_size", h_size)
        if h_size <= 0: w("expected positive integer")
        p("base_bits", h_basebits)
        if h_basebits != 2: w("expected 2 bits per base")
        p("qscore_bits", h_qbits)
        if h_qbits not in (2, 6): w("expected 2 or 6 bits per quality score")
        p("quality_bins", h_bins)

        # Now we have quality binning info, in pairs of 4-byte values.
        # We can also snag the number of tile records which is the next 4 bytes
        qbin_bytes = fh.read((h_bins * 4 * 2) + 4)
        qbin_values = struct.unpack('<'+('II'*h_bins), qbin_bytes[:-4])

        if h_bins:
            for n in range(h_bins):
                print("    bin {n} maps {f}[{f:0{fs}b}] ==> {t}".format(n=n, fs=h_qbits, f=qbin_values[n*2], t=qbin_values[n*2+1]))
        else:
            w("expected to see some bins")

        # Convert bins to map (this works even with no bins)
        qual_bin_map = [0] * (max(qbin_values[::2]) + 1)
        for n in range(h_bins):
            qual_bin_map[qbin_values[n*2]] = qbin_values[n*2+1]

        # Now the number of tile records. We expect 352 for the novaseq
        tile_count, = struct.unpack('<I', qbin_bytes[-4:])
        p("tile_count", tile_count)
        if tile_count != 352: w("expected to see 352 tiles in a novaseq cbcl file")

        # Now the file offsets. I'll unpack these in one go. I think we have 16 bytes per record
        # but the part about the 'non-PF clusters excluded flag' is ambiguous. Is it an extra byte
        # or what?? Hopefully it will become clear.
        all_offset_bytes = fh.read( tile_count * 16 + 1 )
        total_offset = 0
        offset_dict = dict()
        for t in range(tile_count):
            t_number, t_clusters, t_usize, t_csize = struct.unpack('<IIII', all_offset_bytes[t*16:(t+1)*16])

            print("    {t:-3} tile {n} with {c} clusters us={us} cs={cs} off={off}".format(
                        t=t, n=t_number, c=t_clusters, us=t_usize, cs=t_csize, off=total_offset ))

            # Record this vital info so we can access the data for any tile
            offset_dict[t_number] = (total_offset + h_size, t_usize, t_clusters)

            # Tot up the offsets to see where to seek for the next block.
            total_offset += t_csize

        # The documentation says there's another flag:
        #  "non-PF clusters excluded flag -- 1: non-PF clusters are excluded"
        # Oh yes, it's the final byte after the offsets table.

        # Now this is very important because the first 25 cycles are recorded in full (all tiles have
        # 4091904 clusters) but after this the flag goes on and all of the invalid reads are filtered out.
        # So I need to be able to read both types of file, one of which can only be understood with
        # reference to the filter file. How annoying.
        excluded_flag = all_offset_bytes[-1]
        p("excluded_flag", excluded_flag)
        if excluded_flag not in (0, 1): w("excluded_flag can only be 0 or 1")
        excluded_flag = bool(excluded_flag)

        # Now the final total_offset should be the file size plus the header size.
        total_header_size = len(header_bytes) + len(qbin_bytes) + len(all_offset_bytes)

        # And sanity-check agains the claimed header size
        if total_header_size != h_size:
            w("Header claims to be {} bytes but it's actually {}.".format(h_size, total_header_size))

        if total_header_size + total_offset == os.stat(cbcl_file).st_size:
            print("Total file size is {} as expected.".format(total_header_size + total_offset))
        else:
            w("File should be {} bytes but it's actually {}. Difference={}".format(
                        total_header_size + total_offset,
                            os.stat(cbcl_file).st_size,
                                os.stat(cbcl_file).st_size - (total_header_size + total_offset) ))

        if not h_qbits == 2:
            print("The code can't currently read bases where h_qbits != 2")
            exit(0)

        # OK, now let's sample the first and last 100 bases in the first and last 12 tiles. Initially I'll do this
        # without reference to the filter file. Then I'll go add code that makes use of the filter to
        # support excluded_flag=1.
        basemap = ['A', 'C', 'G', 'T']

        for tilenum, (offset, usize, wellcount) in sorted(offset_dict.items())[:12] + sorted(offset_dict.items())[-12:]:

            # Go to the start of the block of tile data and slurp it all
            fh.seek(offset)
            zipdata = gzip.GzipFile(fileobj=fh, mode='rb').read(usize)

            #Print bases in the first 100 and final 100 wells.
            #But I will assume that h_qbits is 2 since apparently it always is (this is now checked above).
            bases_from_start = 100 ; bases_from_end = 100
            seq = []

            # Obtain an indexing array for sequences in a BCL block where the excluded flag is set.
            # Since I'm currently reading a single cycle I need to do this per tile, but in the real
            # code I'll load the offsets once then scan through the .cbcl files per cycle, reading just
            # that tile each time.
            # Obtain this array even for non-excluded .cbcl file because then I can show excluded bases
            # as '-' whereas no-calls within non-excluded sequences will be N's. Got that??
            excluded_offsets, passing_wells = locate_and_load_filter_file(cbcl_file, tilenum)

            # File should correspond with tile layout
            if excluded_flag:
                assert passing_wells == wellcount
                # We can only now discover the real well count
                wellcount = len(excluded_offsets)
            else:
                assert len(excluded_offsets) == wellcount

            for welln in chain(range(bases_from_start), range(wellcount - bases_from_end, wellcount)):

                if excluded_flag:
                    # Indirect offset look-up
                    wellidx = excluded_offsets[welln]
                else:
                    # Direct correspondence between the well number and the file offset
                    # Show the excluded sequences as '-' to be consistent.
                    wellidx = welln if excluded_offsets[welln] != -1 else -1
                    # Or to reveal the base calls for all wells I could do this...
                    # wellidx = welln

                # The manual says "For a two bit quality score, this is two clusters per byte where
                # the bottom 4 bits are the first cluster and the higher 4 bits are the second cluster."
                # TODO - must be 100% sure this is the right way around!! (Can look at the last byte
                # of the odd-numbered filtered data blocks to see that they are like 00001111, as well
                # as comparing my results to the FASTQ files).
                # So the structure of the first byte in a bcl block is:
                #   00 00 00 00 => qual_well_1 call_well_1 qual_well_0 call_well_0
                # By that logic, if any qual is 00 then the call must also be 00 - ie. 00110011 should
                # be impossible, or any other number where there is a pair of zeros not followed by
                # another pair of zeros.
                # I'll add an assertion to sanity-check this.

                if wellidx == -1:
                    # Excluded - though there may still be an underlying basecall
                    # if excluded_flag was False - see above
                    seq.append('-')

                elif wellidx % 2:
                    # Take the high bits
                    if zipdata[wellidx//2] >> 4:
                        assert (zipdata[wellidx//2] >> 6), \
                            "Found base call in high bits with 0 quality: {:08b}".format(zipdata[wellidx//2])
                        seq.append(basemap[ (zipdata[wellidx//2] >> 4) % 4 ])
                    else:
                        seq.append('N')
                else:
                    if zipdata[wellidx//2] % 16:
                        assert (zipdata[wellidx//2] >> 2) % 4, \
                            "Found base call in low bits with 0 quality: {:08b}".format(zipdata[wellidx//2])
                        seq.append(basemap[ zipdata[wellidx//2] % 4 ])
                    else:
                        seq.append('N')

            print( "Sequence in first {} wells of tile {} is {}".format(bases_from_start, tilenum, ''.join(seq[:bases_from_start])) )
            print( "  ...and the last {} wells of tile {} is {}".format(bases_from_end, tilenum, ''.join(seq[bases_from_start:])) )

def locate_and_load_filter_file(cbcl_file, tilenum):
        """ Locate and parse the filter file for this tile.
            Translate this into a list of offsets where the wells will be found in the excluded
            bcl blocks, so if the filter starts 000110101 then the lookup needs to be
            [ -1, -1, -1, 0, 1, -1, 2, -1, 3 ]
            Also return the number of passing wells (final offset + 1)
        """
        # Filter files live one directory up and they have names like s_2_2488.filter ==> s_{lane}_{tile}.filter
        lane_dir = os.path.dirname(os.path.dirname(os.path.realpath(cbcl_file)))

        # lane is 1,2,3 or 4
        lane = lane_dir[-1]

        filter_file = "{d}/s_{lane}_{tile}.filter".format(d=lane_dir, lane=lane, tile=tilenum)

        with open(filter_file, 'rb') as filt_fh:

            filt_header = filt_fh.read(12)
            assert tuple(struct.unpack('<III', filt_header))[:2] == (0, 3)
            num_clusters = struct.unpack('<III', filt_header)[2]

            # Slurp the whole thing - file length should match what the header says.
            filt_bytes = struct.unpack('<{}B'.format(num_clusters), filt_fh.read())

            # Make the table off offsets
            filt_offsets = [-1] * num_clusters
            offset = 0
            for n, flag in enumerate(filt_bytes):
                if flag & 0b00000001:
                    filt_offsets[n] = offset
                    offset += 1

        return filt_offsets, offset

if __name__ == '__main__':
    main(*sys.argv[1:])
