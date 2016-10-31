
#include "fusion.hpp"
#include "process.hpp"
#include "revcomp.hpp"

Fusion::Fusion( Contiglist *contigs, Readlist *reads ) : contigs(contigs), reads(reads){
  fusions_completed = 0;
}

// creates id of fused contigs
std::string Fusion::get_fused_id( std::string contig1_id, std::string contig2_id ){
  if( contig1_id.length() >= 5 && contig1_id.compare( 0, 5, "fused" ) == 0 ){
    contig1_id = contig1_id.substr( 5 );
  }

  if( contig1_id.length() >= 5 && contig2_id.compare( 0, 5, "fused" ) == 0 ){
    contig2_id = contig2_id.substr( 5 );
  }

  return contig1_id+"_<>_"+contig2_id;
}

// complete contig_fusion process
void Fusion::contig_fusion_log( Mismatch fus ){
  // print messages to logfile about current actions
  Log::Inst()->log_it( "Contig fused: " );
  Log::Inst()->log_it( "  Overlap length: " + std::to_string(fus.get_length()) );
  Log::Inst()->log_it( "  Mismatch_score: " + std::to_string(fus.get_score()) );
  Log::Inst()->log_it( "  Contig_i: " + contigs->get_contig(fus.get_index_i()).get_contig_id() );
  Log::Inst()->log_it( "  Contig_j: " + contigs->get_contig(fus.get_index_j()).get_contig_id() );
  Log::Inst()->log_it( "" );
}

// complete contig_fusion process
void Fusion::commit_fusion( std::string fused, std::string fused_id, int index_i, int index_j ){
  contigs->append_contig( 1, contigs->get_contig(index_i) );
  Log::Inst()->log_it( "Contig moved to fused file: " + contigs->get_contig(index_i).get_contig_id() );

  contigs->append_contig( 1, contigs->get_contig(index_j) );
  Log::Inst()->log_it( "Contig moved to fused file: " + contigs->get_contig(index_j).get_contig_id() );

  Log::Inst()->log_it( "Committing: " + fused_id );
  Log::Inst()->log_it( "\t" + fused );
  contigs->append_contig( 0, Contig( reads, fused, fused_id ));

  fusions_completed++;
}

// check overlap section for mismatches
Mismatch Fusion::overlap_check( std::string contig_a, std::string contig_b, int overlap, int end_i, int end_j ){
  Mismatch mim;
  mim.set_end_i( end_i );
  mim.set_end_j( end_j );
  double score_f = 1.0;
  double score_r = 1.0;
  double score = 1.0;

  for( int i=overlap-1; i>=min_overlap; i-- ){
    score_f = mismatch_score( contig_a.substr( contig_a.length() - i, i/2 ), contig_b.substr( 0, i/2 ) );
    score_r = mismatch_score( contig_a.substr( contig_a.length() - (i-i/2) ), contig_b.substr( i/2, i-i/2 ) );
    score = (score_f + score_r)/2;

    // check if this overlap is the best so far
    if( score < mim.get_score() ){
      mim.set_score( score );
      mim.set_length( i );
    }
    // if the score for the end of contig_a is better than the threshold, check the end of contig_b to see if there could be support in the reads for that end going in a different direction
    else if( score_r <= mismatch_threshold ){
      score_f = check_fusion_support( contig_b.substr(i/2, i-i/2), contig_a.substr( contig_a.length() - i, i/2 ) );
      // make sure the support score is at least as good as the threshold
      if( score_f <= mismatch_threshold ){
        score = (score_f + score_r)/2;
        if( score < mim.get_score() ){
          mim.set_score( score );
          mim.set_length( i );
        }
      }
    }
    // if the score for the end of contig_b is better than the threshold, check the end of contig_a to see if there could be support in the reads for that end going in a different direction
    else if( score_f <= mismatch_threshold ){
      score_r = check_fusion_support( contig_a.substr( contig_a.length() - (i-i/2) ), contig_b.substr( i/2, i-i/2 ) );
      // make sure the support score is at least as good as the threshold
      if( score_r <= mismatch_threshold ){
        score = (score_r + score_f)/2;
        if( score < mim.get_score() ){
          mim.set_score( score );
          mim.set_length( i );
        }
      }
    }
  }
  return mim;
}

// create fused contig string
std::string Fusion::build_fusion_string( std::string contig_a, std::string contig_b, int overlap ){
  std::string fused( contig_a.substr( 0, contig_a.length() - (overlap/2) ) );
  fused.append( contig_b.substr( overlap/2 + overlap%2 ) );

  return fused;
}

// remove duplicates from contig remove list
void Fusion::dedup_removals(){
  for( int i=0; i<(int)contig_remove_list.size()-1; i++ ){
    for( int j=i+1; j<contig_remove_list.size(); j++ ){
      if( contig_remove_list[i] == contig_remove_list[j] ){
        contig_remove_list.erase( contig_remove_list.begin() + j );
        j--;
      }
    }
  }
}

// sort index list for removing contigs
void Fusion::sort_removals(){
  int plc_hld;
  for( int i=0; i<(int)contig_remove_list.size()-1; i++ ){
    int low_idx = i;
    for( int j=i+1; j<(int)contig_remove_list.size(); j++ ){
      if( contig_remove_list[j] < contig_remove_list[low_idx] ){
        low_idx = j;
      }
    }
    if( low_idx != i ){
      plc_hld = contig_remove_list[i];
      contig_remove_list[i] = contig_remove_list[low_idx];
      contig_remove_list[low_idx] = plc_hld;
    }
  }
}

// remove fused contigs from contigs list
void Fusion::process_removals(){
  dedup_removals();
  sort_removals();
  for( int i=(int)contig_remove_list.size()-1; i>=0; i-- ){
    contigs->remove_contig( contig_remove_list[i] );
  }
  contig_remove_list.clear();
}

// compile list of best mismatch scores between contigs that meet the mismatch threshold
std::vector<Mismatch> Fusion::get_mismatch_scores( bool first_run ){
  std::vector<Mismatch> match_list;
  int overlap_len = extend_len * 4;

  // loop through each contig to get the end of the contig
  for( int i=0; i<contigs->get_list_size(); i++ ){
    for( int j=i+1; j<contigs->get_list_size(); j++ ){
      std::string contig_i( contigs->get_contig(i).get_sequence() );
      std::string contig_j( contigs->get_contig(j).get_sequence() );
      std::string contig_j_rev( revcomp( contig_j ) );
      int overlap = overlap_len;
      Mismatch mim;

      if( first_run ){
        overlap = (contig_i.length() < contig_j.length()) ? contig_i.length() : contig_j.length();
      }
      else{
        if( overlap > contig_i.length() ){
          overlap = contig_i.length();
        }
        if( overlap > contig_j.length() ){
          overlap = contig_j.length();
        }
      }

      //// Processing of rear end of the contig
      // orientation: i to j
      mim = overlap_check( contig_i, contig_j, overlap, 1, 0 );

      if( mim.get_score() <= mismatch_threshold ){
        mim.set_indices( i, j );
        match_list.push_back( mim );
      }

      // orientation: i to j_rev
      mim = overlap_check( contig_i, contig_j_rev, overlap, 1, 1 );

      if( mim.get_score() <= mismatch_threshold ){
        mim.set_indices( i, j );
        match_list.push_back( mim );
      }

      //// Processing of front end of the contig
      // orientation: j to i
      mim = overlap_check( contig_j, contig_i, overlap, 0, 1);

      if( mim.get_score() <= mismatch_threshold ){
        mim.set_indices( i, j );
        match_list.push_back( mim );
      }

      // orientation: j_rev to i
      mim = overlap_check( contig_j_rev, contig_i, overlap, 0, 0 );

      if( mim.get_score() <= mismatch_threshold ){
        mim.set_indices( i, j );
        match_list.push_back( mim );
      }
    }
  }

  return match_list;
}

// sort the match_list for easier
void Fusion::sort_matches(){
  Mismatch temp_mim;

  // sort by overlap length
  for( int i=0; i<match_list.size(); i++ ){
    for( int j=i; j<match_list.size(); j++ ){
      if( match_list[i].get_length() < match_list[j].get_length() ){
        temp_mim = match_list[i];
        match_list[i] = match_list[j];
        match_list[j] = temp_mim;
      }
    }
  }

  // sort by score
  for( int i=0; i<match_list.size(); i++ ){
    for( int j=i; j<match_list.size(); j++ ){
      if( match_list[i].get_score() > match_list[j].get_score() ){
        temp_mim = match_list[i];
        match_list[i] = match_list[j];
        match_list[j] = temp_mim;
      }
    }
  }
}

// cleans match_list from conflicting matches
void Fusion::clean_matches(){
  for( int i=0; i<(int)match_list.size()-1; i++ ){
    int index_i = match_list[i].get_index_i();
    int end_i = match_list[i].get_end_i();
    int index_j = match_list[i].get_index_j();
    int end_j = match_list[i].get_end_j();
    for( int j=i+1; j<match_list.size(); j++ ){
      if( match_list[j].get_index_i() == index_i && match_list[j].get_end_i() == end_i ){
        match_list.erase( match_list.begin() + j );
        j--;
      }
      else if( match_list[j].get_index_j() == index_i && match_list[j].get_end_j() == end_i ){
        match_list.erase( match_list.begin() + j );
        j--;
      }
      else if( match_list[j].get_index_j() == index_j && match_list[j].get_end_j() == end_j ){
        match_list.erase( match_list.begin() + j );
        j--;
      }
      else if( match_list[j].get_index_i() == index_j && match_list[j].get_end_i() == end_j ){
        match_list.erase( match_list.begin() + j );
        j--;
      }
    }
  }
}

// process the compiled list of fusions
void Fusion::process_fusions(){
  for( int i=0; i<match_list.size(); i++ ){
    // fuse contigs
    int index_i = match_list[i].get_index_i();
    int index_j = match_list[i].get_index_j();
    std::string contig_i = contigs->get_contig(index_i).get_sequence();
    std::string contig_j = contigs->get_contig(index_j).get_sequence();
    std::string fused_id("");
    std::string fused("");
    std::string rev_i("");
    std::string rev_j("");

    // find the reverse compliment of the contigs where necessary
    if( match_list[i].get_end_i() == 0 ){
      contig_i = revcomp( contig_i );
      rev_i = "_r";
    }

    if( match_list[i].get_end_j() == 1 ){
      contig_j = revcomp( contig_j );
      rev_j = "_r";
    }

    contig_fusion_log( match_list[i] );
    fused_id = get_fused_id( contigs->get_contig(index_i).get_contig_id()+rev_i, contigs->get_contig(index_j).get_contig_id()+rev_j );
    fused_id = "fused(" + fused_id + ")";

    // push each index onto the remove vector
    contig_remove_list.push_back( index_i );
    contig_remove_list.push_back( index_j );

    // create fused contig
    fused = build_fusion_string( contig_i, contig_j, match_list[i].get_length() );

    // commit fusion here
    commit_fusion( fused, fused_id, index_i, index_j );
    int new_index = contigs->get_list_size()-1;

    // check for additional appearances of the current contigs in match_list
    for( int j=i+1; j<match_list.size(); j++ ){
      // change reference index and end variables where necessary
      if( index_i == match_list[j].get_index_i() ){
        match_list[j].set_index_i( new_index );
        match_list[j].set_end_i( 0 );
      }
      else if( index_i == match_list[j].get_index_j() ){
        match_list[j].set_index_j( new_index );
        match_list[j].set_end_j( 0 );
      }
      else if( index_j == match_list[j].get_index_i() ){
        match_list[j].set_index_i( new_index );
        match_list[j].set_end_i( 1 );
      }
      else if( index_j == match_list[j].get_index_j() ){
        match_list[j].set_index_j( new_index );
        match_list[j].set_end_j( 1 );
      }
    }
  }
}

// contig_fusion: Attempt to support fusion in case of possibly poorly constructed end.. returns new score from section in question
//    ::> contig object is the second while contig_ref is the first and the extension is being made off the front of the object
double Fusion::check_fusion_support( std::string contig, std::string contig_ref ){
  Extension *exten = new Extension( reads, contig_ref.length(), contig );

  std::string support_string( "" );
  double score = 0;
  const int pos = contig_ref.length() - 1;
  const int start = -1;

  exten->matches.start_match();

  // return if matches found is less than min_cov
  if( exten->matches.get_matchlist_size() < min_cov ){
    delete exten;
    return 1.0;
  }

  // create missed bp's vector to keep track of how many bp's each read contains that are below the max at that position
  std::vector<int> mismatch( exten->matches.get_matchlist_size(), 0 );
  std::vector<int> ambiguous_bp( exten->matches.get_matchlist_size(), 0 );
  int mismatch_tot = 0;
  int mismatch_avg = 0;
  int cmp_len = 0;

  // count mismatch's in each read match
  for( int i=0; i<contig_ref.length(); i++ ){
    int next_char_ref = contig_ref[pos-i];
    for( int j=0; j<exten->matches.get_matchlist_size(); j++ ){
      int next_char = exten->matches.get_pos( j, start-i );
      if( next_char == 'N' ){
        ambiguous_bp[j]++;
      }
      // check if next_char matches or exists in the current read
      else if( next_char != next_char_ref && next_char != -1 ){
        mismatch[j]++;
      }
    }
  }

  // remove reads with more than the max_missed
  exten->set_missed_bp( mismatch );
  exten->set_missed_avg( max_missed );
  exten->error_removal();

  // remove reads with more than 3 N's
  exten->set_missed_bp( ambiguous_bp );
  exten->set_missed_avg( 3 );
  exten->error_removal();

  // if matchlist is smaller than min_cov, bail, return 1.0 (the least supportive return value)
  if( exten->matches.get_matchlist_size() < min_cov ){
    delete exten;
    return 1.0;
  }

  // get string built from remaining matches.
  exten->build_string();
  support_string = exten->get_extension();
  cmp_len = support_string.length();

  // use the smaller of contig_ref.length() and min_overlap/2 as min_length to compare support_string with matching contig
  if( cmp_len < contig_ref.length() && cmp_len < (min_overlap/2 + min_overlap%2) ){
    delete exten;
    return 1.0;
  }

  // get smallest of string lengths
  if( cmp_len > contig_ref.length() ){
    cmp_len = contig_ref.length();
  }

  score = mismatch_score( support_string.substr( support_string.length() - cmp_len), contig_ref.substr( contig_ref.length() - cmp_len ) );

  delete exten;
  return score;
}

// fuse contigs wherever possible
void Fusion::run_fusion( bool first_run ){
  Log::Inst()->log_it( "Fuse contigs" );

  // MISMATCH SCORES //
  match_list = get_mismatch_scores( first_run );

  // SORT AND CLEAN MATCH LIST //
  sort_matches();
  clean_matches();

  if( verbose ){
    for( int i=0; i<match_list.size(); i++ ){
      Log::Inst()->log_it( "score: " + std::to_string( match_list[i].get_score() ) );
      Log::Inst()->log_it( "i: " + std::to_string( match_list[i].get_index_i() ) +  " j: " + std::to_string( match_list[i].get_index_j() ) );
    }
  }

  // FUSE CONTIGS //
  process_fusions();

  // remove fused contigs from the vector
  process_removals();
}

// tally mismatches in substrings passed and return score in the form of misatches per length
double Fusion::mismatch_score( std::string contig_sub_a, std::string contig_sub_b ){
  int mismatch = 0;

  // protect against division by zero
  if( contig_sub_a.length() == 0 ){
    return 1.0;
  }

  for( int i=0; i<contig_sub_a.length(); i++ ){
    if( contig_sub_a[i] != contig_sub_b[i] ){
      mismatch++;
    }
    // bail when too many mismatches have been tallied.. avoid excessive processing
    if( double(mismatch) /contig_sub_a.length() > mismatch_threshold ){
      return 1.0;
    }
  }

  double score = double(mismatch) / contig_sub_a.length();
  return score;
}
