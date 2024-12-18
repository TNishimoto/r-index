#pragma once
// #include "./dynamic_rlbwt_helper.hpp"
#include "../module/stool/include/stool.hpp"

namespace stool
{
    namespace r_index
    {
        ////////////////////////////////////////////////////////////////////////////////
        /// @class      SamplingSATBuilder
        /// @brief      A builder of sampling suffix arrays
        ///
        ////////////////////////////////////////////////////////////////////////////////
        class SamplingSATBuilder
        {
        public:
            static tuple<string, vector<pair<uint64_t, uint64_t>>, vector<uint64_t>> build(std::string file_path, int message_paragraph = stool::Message::SHOW_MESSAGE)
            {

                using RLBWT = stool::rlbwt2::RLE<uint8_t>;
                // stool::rlbwt2::BWTAnalysisResult analyzer;
                // stool::rlbwt2::RLE<uint8_t> static_rlbwt = stool::rlbwt2::RLE<uint8_t>::build(file_path);
                stool::rlbwt2::RLE<uint8_t> static_rlbwt = stool::rlbwt2::RLE<uint8_t>::build_from_file(file_path, stool::Message::SHOW_MESSAGE);

                stool::WT wt = stool::rlbwt2::WaveletTreeOnHeadChars::build(&static_rlbwt);
                stool::rlbwt2::LightFPosDataStructure fpos_array;
                fpos_array.build(static_rlbwt.get_head_char_vec(), *static_rlbwt.get_lpos_vec(), &wt, stool::Message::SHOW_MESSAGE);
                using LF_DATA = stool::rlbwt2::LFDataStructureBasedOnRLBWT<RLBWT, stool::rlbwt2::LightFPosDataStructure>;
                LF_DATA rle_wt(&static_rlbwt, &fpos_array);

                uint64_t end_marker_lposition = static_rlbwt.get_end_rle_lposition();
                uint64_t end_marker_position = static_rlbwt.get_lpos(end_marker_lposition);

                vector<pair<uint64_t, uint64_t>> samples_first; // text positions corresponding to first characters in BWT runs, and their ranks 0...R-1
                vector<uint64_t> samples_last;                  // text positions corresponding to last characters in BWT runs
                std::string bwt_s;
                bwt_s.resize(static_rlbwt.str_size());

                samples_first.resize(static_rlbwt.rle_size());
                samples_last.resize(static_rlbwt.rle_size());

                // using LF_DATA = stool::rlbwt2::LFDataStructure<stool::rlbwt2::RLE<uint8_t>, stool::rlbwt2::LightFPosDataStructure>;

                // LF_DATA rle_wt(&static_rlbwt, &fpos_array);
                stool::bwt::BackwardISA<LF_DATA> isa_ds;

                // uint64_t p3 = rle_wt.lf(end_marker_lposition);

                uint64_t text_size = static_rlbwt.str_size();

                isa_ds.set(&rle_wt, end_marker_lposition, text_size);
                int64_t text_position = text_size;
                uint8_t min_c = static_rlbwt.get_smallest_character();

                uint64_t message_counter = 10000000;
                uint64_t processed_text_length = 0;
                auto _end = isa_ds.end();

                for (stool::bwt::BackwardISA<LF_DATA>::iterator it = isa_ds.begin(); it != _end; ++it)
                {

                    text_position--;

                    message_counter++;
                    processed_text_length++;

                    if (message_paragraph >= 0 && message_counter > 10000000)
                    {

                        std::cout << stool::Message::get_paragraph_string(message_paragraph + 1) << "Processing... [" << (processed_text_length / 1000000) << "/" << (text_size / 1000000) << "MB] \r" << std::flush;
                        message_counter = 0;
                    }

                    uint64_t lindex = static_rlbwt.get_lindex_containing_the_position(*it);
                    uint64_t run_length = static_rlbwt.get_run(lindex);
                    uint64_t starting_position = static_rlbwt.get_lpos(lindex);
                    uint64_t diff = (*it) - starting_position;
                    uint8_t c = static_rlbwt.get_char_by_run_index(lindex);

                    /*
                    uint64_t lindex = static_rlbwt.get_lindex_containing_the_position(*it);
                    uint64_t run_length = static_rlbwt.get_run(lindex);
                    uint64_t starting_position = static_rlbwt.get_lpos(lindex);
                    uint8_t c = static_rlbwt.get_char_by_run_index(lindex);
                    uint64_t diff = (*it) - starting_position;
                    */

                    bwt_s[*it] = c == min_c ? 1 : c;

                    if (diff == 0)
                    {
                        samples_first[lindex] = std::pair<uint64_t, uint64_t>(text_position, lindex);
                    }
                    if (diff + 1 == run_length)
                    {
                        samples_last[lindex] = text_position;
                    }
                }
                if (message_paragraph >= 0 && text_size > 0)
                {
                    std::cout << std::endl;
                    std::cout << "[END]" << std::endl;
                }

                return std::tuple<string, vector<pair<uint64_t, uint64_t>>, vector<uint64_t>>(bwt_s, samples_first, samples_last);
            }
        };
    }
}
