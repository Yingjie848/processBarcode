#include <iostream>
#include <string> /* for string */
#include <fstream> /* for file */
#include <stdlib.h> /* for exit() */
#include <istream> /* for getopt */
#include <sstream> /* for stream */
#include <vector>
#include <map>
#include <boost/lexical_cast.hpp> /* convert nucmeric to string */
#include "ssw_cpp.h"

#define DEBUG 

using namespace std;

// Declare structure
struct ANNOTATE_T{
    short int type; //1 or 2
    int count_P;
    int count_D;
};

/* Predefined barcode and sub-barcode sequences */
// Please modify these information if needed
string pBarcode = "CCCTAGCGTAACTCTCGAGGTAGTA";
string dBarcode = "CCCTAGCGTAACTCTCGAGACGACG";
int sub_barcode_start=14; //xhoI + barcode 1 (or 2)
int sub_barcode_length=11;
string sub_pBarcode=pBarcode.substr(sub_barcode_start,sub_barcode_length);
string sub_dBarcode=dBarcode.substr(sub_barcode_start,sub_barcode_length);
string xhoI="TCGAG";
string telomere1="TTAGGG";
string telomere2="TAGGGT";
string telomere3="AGGGTT";
string telomere4="GGGTTA";
string telomere5="GGTTAG";
string telomere6="GTTAGG";

/* Declares functions */
static void usage();
map<string,int> readExceptionFile(string exception_list);
static void PrintAlignment(const StripedSmithWaterman::Alignment& alignment);
static void PrintAlignmentTable(const StripedSmithWaterman::Alignment& alignment);
static void PrintAlignmentTable(const StripedSmithWaterman::Alignment& alignment, ofstream &stat);
StripedSmithWaterman::Alignment align_query_to_reference(string barcode, string sub_barcode, string reference);
ANNOTATE_T annotate_seq(int cutoff_match, int cutoff_mismatch, int cutoff_soft_clipping, int matched_bps_1, int matched_bps_2, float matched_perc_1, float matched_perc_2, StripedSmithWaterman::Alignment alignment1, StripedSmithWaterman::Alignment alignment2);
bool is_file_exist(string fileName);

/* main function */
int main( int argc, char** argv ){

    //input parameters
    if( argc < 2 ){
        usage();
    }
    
    int cutoff_match=9;
    int cutoff_mismatch=1;
    int cutoff_soft_clipping=1;
    string seq1_exception_list; 
    string seq2_exception_list; 
    string output_prefix="out";

    int c;
    while ((c=getopt(argc,argv,"m:M:s:a:b:p:")) != -1)
    {
        switch (c)
        {
            case 'm' : if(optarg) cutoff_match = atoi(optarg); break; 
            case 'M' : if(optarg) cutoff_mismatch = atoi(optarg); break; 
            case 's' : if(optarg) cutoff_soft_clipping = atoi(optarg); break; 
            case 'a' : if(optarg) seq1_exception_list=optarg; break;
            case 'b' : if(optarg) seq2_exception_list=optarg; break;
            case 'p' : if(optarg) output_prefix=optarg; break;
        }
    }
    string seq1_file = argv[optind++];
    string seq2_file = argv[optind++];

    string out_seq1_mapped_reads = output_prefix+"_R1.map.fa";
    string out_seq2_mapped_reads = output_prefix+"_R2.map.fa";
    string out_seq1_xhoI_reads = output_prefix+"_R1.map.xhoI.fa";
    string out_seq2_xhoI_reads = output_prefix+"_R2.map.xhoI.fa";
    string out_seq1_telo_reads = output_prefix+"_R1.map.telo.fa";
    string out_seq2_telo_reads = output_prefix+"_R2.map.telo.fa";
    string out_PD1 = output_prefix+".unmap.PD1.fa"; //Output proximal and distal sequence with one end. Only P-D pattern is output. For D-P pattern, sequence will be reversed and complemented. P-P pattern has short insert size will not be output.
    string out_PD2 = output_prefix+".unmap.PD2.fa"; //Output proximal and distal sequence with another end. Only P-D pattern is output. For D-P pattern, sequence will be reversed and complemented. P-P pattern has short insert size will not be output.
    string out_PP1 = output_prefix+".unmap.PP1.fa"; //Output proximal and proximal sequence with one end. 
    string out_PP2 = output_prefix+".unmap.PP2.fa"; //Output proximal and proximal sequence with another end. 
    string out_SPB = output_prefix+".unmap.SPB.fa"; //Output Sequence labeled with Proximal Barcode (SPB)
    string out_SDB = output_prefix+".unmap.SDB.fa"; //Output Sequence labeled with Distal Barcode (SDB)
    string out_stat = output_prefix+".stat";

    //output stat file
    ofstream stat ( out_stat.c_str() );
    if ( ! stat )
        cerr << "fail to create stat file" << out_stat <<endl;

    stat << "#Parameters: " << " -m " << cutoff_match << " -M " << cutoff_mismatch 
         << " -a " << seq1_exception_list << " -b " << seq2_exception_list << " -p " << output_prefix 
         << " R1.fa: " << seq1_file << " R2.fa: " << seq2_file; 
    
    // Declare
    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment1;
    StripedSmithWaterman::Alignment alignment2;
    StripedSmithWaterman::Alignment alignment3;
    StripedSmithWaterman::Alignment alignment4;
    // Declares map
    map<string,int> map_seq1_exception;
    map<string,int> map_seq2_exception;
    map<string,int> map_pBarcode; 
    map<string,int>::iterator map_pBarcode_Iter; 
    map<string,int> map_dBarcode; 
    map<string,int>::iterator map_dBarcode_Iter; 
    map<string,int> map_pattern; 
    map<string,int>::iterator map_pattern_Iter; 
    //
    string seq1;  //input sequence
    string seq2;  //input sequence
    int count_seq=0;       //count total # of reads
    int count_seq1_exception=0; //count # of excepted reads
    int count_seq2_exception=0; //count # of excepted reads
    int count_seq1_xhoI=0;      //count xhoI
    int count_seq1_telo=0;      //count telomere
    int count_seq2_xhoI=0;      //count xhoI
    int count_seq2_telo=0;      //count telomere
    int count_seq_unmapped=0; //count unmapped reads
    int count_seq1_P=0;    //count R1 proximal reads
    int count_seq1_D=0;    //count R1 distal reads
    int count_seq2_P=0;    //count R2 proximal reads
    int count_seq2_D=0;    //count R2 distal reads
    int count_seq1_usable_proximal_barcode=0;  //count usable proximal barcode, include P-P, P-D, D-P,P-None, None-P
    int count_seq2_usable_proximal_barcode=0;  //count usable proximal barcode, include P-P, P-D, D-P,P-None, None-P
    int count_seq1_usable_distal_barcode=0;  //count usable distal barcode, include P-D, D-P
    int count_seq2_usable_distal_barcode=0;  //count usable distal barcode, include P-D, D-P

    //open files
    ifstream infile1;
    infile1.open(seq1_file.c_str());
    if ( ! infile1 )
        cerr << "fail to open file" << seq1_file << endl;
     
    ifstream infile2;
    infile2.open(seq2_file.c_str());
    if ( ! infile2 )
        cerr << "fail to open file" << seq2_file << endl;

    //output files
    ofstream map_seq1 ( out_seq1_mapped_reads.c_str() );
    if ( ! map_seq1 )
        cerr << "fail to create read1 clean file" << out_seq1_mapped_reads <<endl;

    ofstream map_seq2 ( out_seq2_mapped_reads.c_str() );
    if ( ! map_seq2 )
        cerr << "fail to create read1 clean file" << out_seq2_mapped_reads <<endl;

    ofstream map_seq1_xhoI ( out_seq1_xhoI_reads.c_str() );
    if ( ! map_seq1_xhoI )
        cerr << "fail to create read1 clean file" << out_seq1_xhoI_reads <<endl;

    ofstream map_seq2_xhoI ( out_seq2_xhoI_reads.c_str() );
    if ( ! map_seq2_xhoI )
        cerr << "fail to create read1 clean file" << out_seq2_xhoI_reads <<endl;

    ofstream map_seq1_telo ( out_seq1_telo_reads.c_str() );
    if ( ! map_seq1_telo )
        cerr << "fail to create read1 clean file" << out_seq1_telo_reads <<endl;

    ofstream map_seq2_telo ( out_seq2_telo_reads.c_str() );
    if ( ! map_seq2_telo )
        cerr << "fail to create read1 clean file" << out_seq2_telo_reads <<endl;

    ofstream clean1 ( out_PD1.c_str() );
    if ( ! clean1 )
        cerr << "fail to create read1 clean file" << out_PD1 <<endl;
    
    ofstream clean2 ( out_PD2.c_str() );
    if ( ! clean2 )
        cerr << "fail to create read1 clean file" << out_PD2 <<endl;
    
    ofstream clean3 ( out_PP1.c_str() );
    if ( ! clean3 )
        cerr << "fail to create read1 clean file" << out_PP1 <<endl;
    
    ofstream clean4 ( out_PP2.c_str() );
    if ( ! clean4 )
        cerr << "fail to create read1 clean file" << out_PP2 <<endl;

    ofstream clean5 ( out_SPB.c_str() );
    if ( ! clean5 )
        cerr << "fail to create read1 clean file" << out_SPB <<endl;
    
    ofstream clean6 ( out_SDB.c_str() );
    if ( ! clean6 )
        cerr << "fail to create read1 clean file" << out_SDB <<endl;

    //load exception list
    if( is_file_exist(seq1_exception_list) )
        map_seq1_exception=readExceptionFile(seq1_exception_list);

    if( is_file_exist(seq2_exception_list) )
        map_seq2_exception=readExceptionFile(seq2_exception_list);


    while( getline( infile1, seq1, '\n' ) ) //read line by '\n'
    {
            getline( infile2, seq2, '\n' );
            count_seq++; //count total sequence
            if( map_seq1_exception[boost::lexical_cast<string>(count_seq-1)]==1 || map_seq2_exception[boost::lexical_cast<string>(count_seq-1)]==1 )
            {
                //cout << count_seq << " " <<  map_seq1_exception[boost::lexical_cast<string>(count_seq-1)] << " " <<  map_seq2_exception[boost::lexical_cast<string>(count_seq-1)] << endl;
            if( map_seq1_exception[boost::lexical_cast<string>(count_seq-1)]==1 )
            {
                count_seq1_exception++;
                map_seq1 << seq1 << endl;

                //align xhoI
                if( seq1.find(xhoI) == 0 ){
                    count_seq1_xhoI++;
                    map_seq1_xhoI << seq1 << endl;
                }else if( (seq1.find(telomere1) == 0) || (seq1.find(telomere2) == 0) || (seq1.find(telomere3) == 0) || (seq1.find(telomere4) == 0) || (seq1.find(telomere5) == 0) || (seq1.find(telomere6) == 0) ){
                    count_seq1_telo++;
                    map_seq1_telo << seq1 << endl;
                } 
                continue;
            }
            if( map_seq2_exception[boost::lexical_cast<string>(count_seq-1)]==1 )
            {
                count_seq2_exception++;
                map_seq2 << seq2 << endl;

                //align xhoI
                if( seq2.find(xhoI) == 0 ){
                    count_seq2_xhoI++;
                    map_seq2_xhoI << seq2 << endl;
                }else if( (seq2.find(telomere1) == 0) || (seq2.find(telomere2) == 0) || (seq2.find(telomere3) == 0) || (seq2.find(telomere4) == 0) || (seq2.find(telomere5) == 0) || (seq2.find(telomere6) == 0) ){
                    count_seq2_telo++;
                    map_seq2_telo << seq1 << endl;
                } 
                continue;
            }
            }
#ifdef DEBUG
            cout << "#" << count_seq <<endl; 
#endif            
            count_seq_unmapped++;
            /*     process seq1     */
            //align proximal barcode, search barcode using string::find function, then do alignment using smith-waterman method
            alignment1 = align_query_to_reference(pBarcode,sub_pBarcode,seq1);
            int matched_bps_1=(alignment1.ref_end-alignment1.ref_begin+1)-alignment1.mismatches;
            float matched_perc_1=(float)matched_bps_1*100/(alignment1.ref_end-alignment1.ref_begin+1);
            
            //align distal barcode
            alignment2 = align_query_to_reference(dBarcode,sub_dBarcode,seq1); 
            int matched_bps_2=(alignment2.ref_end-alignment2.ref_begin+1)-alignment2.mismatches;
            float matched_perc_2=(float)matched_bps_2*100/(alignment2.ref_end-alignment2.ref_begin+1);
#ifdef DEBUG
            //cout << seq1 << endl;
            cout << "P:" << matched_bps_1 << "/" << matched_perc_1 << "%\tD:" << matched_bps_2 << "/" << matched_perc_2 << "%" << endl;
#endif

            //annotate
            //filter by matched base pairs, mismatch number, reference start of alignment, query end of alignment, percent of matched base pairs
            ANNOTATE_T annotate_seq1 = annotate_seq(cutoff_match, cutoff_mismatch, cutoff_soft_clipping, matched_bps_1, matched_bps_2, matched_perc_1, matched_perc_2, alignment1, alignment2);
            count_seq1_P+=annotate_seq1.count_P;
            count_seq1_D+=annotate_seq1.count_D;

            /*     process seq2     */
            //align proximal barcode
            alignment3 = align_query_to_reference(pBarcode,sub_pBarcode,seq2); 
            int matched_bps_3=(alignment3.ref_end-alignment3.ref_begin+1)-alignment3.mismatches;
            float matched_perc_3=(float)matched_bps_3*100/(alignment3.ref_end-alignment3.ref_begin+1);
            
            //align distal barcode
            alignment4 = align_query_to_reference(dBarcode,sub_dBarcode,seq2); 
            int matched_bps_4=(alignment4.ref_end-alignment4.ref_begin+1)-alignment4.mismatches;
            float matched_perc_4=(float)matched_bps_4*100/(alignment4.ref_end-alignment4.ref_begin+1);
#ifdef DEBUG
            //cout << seq2 << endl; 
            cout << "P:" << matched_bps_3 << "/" << matched_perc_3 << "%\tD:" << matched_bps_4 << "/" << matched_perc_4 << "%" << endl;
#endif

            //annotate
            ANNOTATE_T annotate_seq2 = annotate_seq(cutoff_match, cutoff_mismatch, cutoff_soft_clipping, matched_bps_3, matched_bps_4, matched_perc_3, matched_perc_4, alignment3, alignment4); 
            count_seq2_P+=annotate_seq2.count_P;
            count_seq2_D+=annotate_seq2.count_D;

#ifdef DEBUG
            //print annotated barcode patterns
            cout << "Pattern:\t" << annotate_seq1.type << "-" << annotate_seq2.type << endl;
#endif
            map_pattern[boost::lexical_cast<string>(annotate_seq1.type)+"-"+boost::lexical_cast<string>(annotate_seq2.type)]++;
            
            //output clean sequence
            if( (annotate_seq1.type == 1 && annotate_seq2.type == 1) ) //P-P
            { 
                map_pBarcode[seq1.substr(alignment1.ref_begin,alignment1.ref_end-alignment1.ref_begin+1)]++;
                clean3 << seq1.substr(alignment1.ref_end+1,seq1.size()-alignment1.ref_end+1) << endl;
                clean5 << seq1.substr(alignment1.ref_end+1,seq1.size()-alignment1.ref_end+1) << endl;
                map_pBarcode[seq2.substr(alignment3.ref_begin,alignment3.ref_end-alignment3.ref_begin+1)]++;
                clean4 << seq2.substr(alignment3.ref_end+1,seq1.size()-alignment3.ref_end+1) << endl;
                clean5 << seq2.substr(alignment3.ref_end+1,seq1.size()-alignment3.ref_end+1) << endl;
                count_seq1_usable_proximal_barcode++;
                count_seq2_usable_proximal_barcode++;
            }
            else if( (annotate_seq1.type == 1 && annotate_seq2.type == 2)   //P-D
                     || (annotate_seq1.type == 2 && annotate_seq2.type == 1) ) //D-P
            {
                if( annotate_seq1.type == 1 )      //clean seq1 P
                {
                    map_pBarcode[seq1.substr(alignment1.ref_begin,alignment1.ref_end-alignment1.ref_begin+1)]++;    
                    clean1 << seq1.substr(alignment1.ref_end+1,seq1.size()-alignment1.ref_end+1) << endl;
                    clean5 << seq1.substr(alignment1.ref_end+1,seq1.size()-alignment1.ref_end+1) << endl;
                    count_seq1_usable_proximal_barcode++;
                }
                else if( annotate_seq1.type == 2 ) //clean seq1 D
                {
                    map_dBarcode[seq1.substr(alignment2.ref_begin,alignment2.ref_end-alignment2.ref_begin+1)]++;    
                    clean2 << seq1.substr(alignment2.ref_end+1,seq1.size()-alignment2.ref_end+1) << endl;
                    clean6 << seq1.substr(alignment2.ref_end+1,seq1.size()-alignment2.ref_end+1) << endl;
                    count_seq1_usable_distal_barcode++;
                }
                if( annotate_seq2.type == 1 )      //clean seq2 P
                {    
                    map_pBarcode[seq2.substr(alignment3.ref_begin,alignment3.ref_end-alignment3.ref_begin+1)]++;
                    clean1 << seq2.substr(alignment3.ref_end+1,seq1.size()-alignment3.ref_end+1) << endl;
                    clean5 << seq2.substr(alignment3.ref_end+1,seq1.size()-alignment3.ref_end+1) << endl;
                    count_seq2_usable_proximal_barcode++;
                }
                else if( annotate_seq2.type == 2 ) //clean seq2 D
                {
                    map_dBarcode[seq2.substr(alignment4.ref_begin,alignment4.ref_end-alignment4.ref_begin+1)]++;
                    clean2 << seq2.substr(alignment4.ref_end+1,seq1.size()-alignment4.ref_end+1) << endl;
                    clean6 << seq2.substr(alignment4.ref_end+1,seq1.size()-alignment4.ref_end+1) << endl;
                    count_seq2_usable_distal_barcode++;
                }
            }
            else if( (annotate_seq1.type == 0 && annotate_seq2.type == 1 )      //None-P
                     || (annotate_seq1.type == 1 && annotate_seq2.type == 0) ) //P-None
            {
                if( annotate_seq1.type == 1 )      //clean seq1 P
                {
                    map_pBarcode[seq1.substr(alignment1.ref_begin,alignment1.ref_end-alignment1.ref_begin+1)]++;    
                    clean5 << seq1.substr(alignment1.ref_end+1,seq1.size()-alignment1.ref_end+1) << endl;
                    count_seq1_usable_proximal_barcode++;
                }
                if( annotate_seq2.type == 1 )      //clean seq2 P
                {
                    map_pBarcode[seq2.substr(alignment3.ref_begin,alignment3.ref_end-alignment3.ref_begin+1)]++;
                    clean5 << seq2.substr(alignment3.ref_end+1,seq2.size()-alignment3.ref_end+1) << endl;
                    count_seq2_usable_proximal_barcode++;
                }
            }
            else if( (annotate_seq1.type == 0 && annotate_seq2.type == 2)     //None-D
                     || (annotate_seq1.type == 2 && annotate_seq2.type == 0)  //D-None
                     || (annotate_seq1.type == 2 && annotate_seq2.type == 2) )//D-D
            {
                /*
                if( annotate_seq1.type == 2 )      //clean seq1 D
                {
                    map_pBarcode[seq1.substr(alignment2.ref_begin,alignment2.ref_end-alignment2.ref_begin+1)]++;    
                    clean6<< seq1.substr(alignment2.ref_end+1,seq1.size()-alignment2.ref_end+1) << endl;
                }
                if( annotate_seq2.type == 2 )      //clean seq2 D
                {
                    map_pBarcode[seq2.substr(alignment4.ref_begin,alignment4.ref_end-alignment4.ref_begin+1)]++;    
                    clean6<< seq2.substr(alignment4.ref_end+1,seq2.size()-alignment4.ref_end+1) << endl;
                }
                */
            }
            else{
                ;
            }
    }
    //Output annotated barcode sequence (include mismatches)
    for( map_pBarcode_Iter=map_pBarcode.begin();map_pBarcode_Iter!=map_pBarcode.end();map_pBarcode_Iter++)
        stat << "Proximal " << map_pBarcode_Iter->first << "\t" << map_pBarcode_Iter->second << endl;
    for( map_dBarcode_Iter=map_dBarcode.begin();map_dBarcode_Iter!=map_dBarcode.end();map_dBarcode_Iter++)
        stat << "Distal   " << map_dBarcode_Iter->first << "\t" << map_dBarcode_Iter->second << endl;
    stat << endl;
    //Output annotated barcode patterns
    stat << "Barcode patterns (0:none, 1:R1 P, 2:R1 D, 3:R2 P, 4:R2 D)\n";
    for( map_pattern_Iter=map_pattern.begin();map_pattern_Iter!=map_pattern.end();map_pattern_Iter++)
        stat << map_pattern_Iter->first << "\t" << map_pattern_Iter->second << endl;
    stat << endl;
    //count annotated sequence
    stat << "Reads\tTotal_reads\tException\tXhoI_reads\tTelomere_reads\tTotal_reads_unmapped\tProximal_barcode_annotated\tDistal_barcode_annotated\tUsable_proximal_barcode\tUsable_distal_barcode\n"; 
    stat << "R1\t" << count_seq << "\t" << count_seq1_exception << "\t" << count_seq1_xhoI << "\t" << count_seq1_telo << "\t" << count_seq_unmapped << "\t" << count_seq1_P << "\t" << count_seq1_D << "\t" << count_seq1_usable_proximal_barcode << "\t" << count_seq1_usable_distal_barcode << endl;
    stat << "R2\t" << count_seq << "\t" << count_seq2_exception << "\t" << count_seq2_xhoI << "\t" << count_seq2_telo << "\t" << count_seq_unmapped << "\t" << count_seq2_P << "\t" << count_seq2_D << "\t" << count_seq2_usable_proximal_barcode << "\t" << count_seq2_usable_distal_barcode << endl;

    infile1.close();infile2.close();
    clean1.close();clean2.close();clean3.close();
    clean4.close();clean5.close();clean6.close();

    return 0;
}

/* ==========USAGE function=========== */
static void usage ()
{	
    cout << "Usage: process_barcode -m 10 -M 1 -s 1 -a <R1_exception_list> -b <R2_exception_list> -p <output_prefix>  <R1.fa> <R2.fa>\n"
         << " -m <int>      minimum of match number, default is 9\n"
         << " -M <int>      max mismatch number, default is 1\n"
         << " -s <int>      max soft clipping base pairs, default is 1\n"
         << " -a <string>,  Exception list for R1\n"
         << " -b <string>,  Exception list for R2\n"
         << " -p <string>,  output prefix for all output files\n"
         << " <R1.fa>       input sequence R1\n"
         << " <R2.fa>       input sequence R2\n";
    exit(1);
}

/* ==========PrintAlignment function========== */
static void PrintAlignment(const StripedSmithWaterman::Alignment& alignment){
  cout << "=== SSW result ===" << endl;
  cout << "Best Smith-Waterman score:\t" << alignment.sw_score << endl
       << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << endl
       << "Reference start:\t" << alignment.ref_begin << endl
       << "Reference end:\t" << alignment.ref_end << endl
       << "Query start:\t" << alignment.query_begin << endl
       << "Query end:\t" << alignment.query_end << endl
       << "Next-best reference end:\t" << alignment.ref_end_next_best << endl
       << "Number of mismatches:\t" << alignment.mismatches << endl
       << "Cigar: " << alignment.cigar_string << endl;
  cout << "===" << endl;
}

/* ==========PrintAlignmentTable function========== */
static void PrintAlignmentTable(const StripedSmithWaterman::Alignment& alignment){
  cout << "<- Aln: " ;
  cout << "  " << alignment.sw_score << "\t" << alignment.ref_begin << "\t" << alignment.ref_end << "\t" 
       << alignment.query_begin << "\t" << alignment.query_end << "\t" << alignment.mismatches << "\t"  
       << alignment.cigar_string;
  cout << " ->" << endl;
}
/* ==========PrintAlignmentTable function========== */
static void PrintAlignmentTable(const StripedSmithWaterman::Alignment& alignment, ofstream &stat){
  stat << "<-Aln: ";
  stat << "  " << alignment.sw_score << "\t" << alignment.ref_begin << "\t" << alignment.ref_end << "\t" 
       << alignment.query_begin << "\t" << alignment.query_end << "\t" << alignment.mismatches << "\t"  
       << alignment.cigar_string << endl;
  stat << "->" << endl;
}

StripedSmithWaterman::Alignment align_query_to_reference(string barcode, string sub_barcode, string reference){
    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;
    
    if( reference.find(barcode) == 0 )
    {   
        alignment.ref_begin=0;
        alignment.ref_end=barcode.size()-1;
        alignment.query_begin=0;
        alignment.query_end=barcode.size()-1;
        alignment.mismatches=0;
        alignment.cigar_string=boost::lexical_cast<string>(barcode.size())+"=";
        alignment.sw_score=0;
     }
     else if( reference.find(sub_barcode) == 0 )
     {   
        alignment.ref_begin=0;
        alignment.ref_end=sub_pBarcode.size()-1;
        alignment.query_begin=sub_barcode_start;
        alignment.query_end=sub_barcode_start+sub_barcode_length-1;
        alignment.mismatches=0;
        alignment.cigar_string=boost::lexical_cast<string>(sub_barcode_length)+"=";
        alignment.sw_score=0;
     }
     else{
        aligner.Align(barcode.c_str(), reference.c_str(), reference.size(), filter, &alignment);
     }
     return alignment;
}


ANNOTATE_T annotate_seq(int cutoff_match, int cutoff_mismatch, int cutoff_soft_clipping, 
     int matched_bps_1, int matched_bps_2, float matched_perc_1, float matched_perc_2, 
     StripedSmithWaterman::Alignment alignment1, StripedSmithWaterman::Alignment alignment2 )
{
            
            ANNOTATE_T annotate;
            annotate.type=0;
            annotate.count_P=0;
            annotate.count_D=0;
            if( matched_bps_1 >= cutoff_match && alignment1.mismatches <= cutoff_mismatch 
                && alignment1.ref_begin <= cutoff_mismatch && alignment1.query_end >= pBarcode.size()-cutoff_soft_clipping-1 
                && matched_perc_1 > matched_perc_2 )
            {
#ifdef DEBUG
                cout << "P: " << pBarcode << endl;
                PrintAlignmentTable(alignment1);
#endif
                annotate.count_P++;
                annotate.type=1;
            }
            else if( matched_bps_2 >= cutoff_match && alignment2.mismatches <= cutoff_mismatch 
                && alignment2.ref_begin <= cutoff_mismatch && alignment2.query_end >= pBarcode.size()-cutoff_soft_clipping-1 
                && matched_perc_1 < matched_perc_2 )
            {
#ifdef DEBUG
                cout << "D: " << dBarcode << endl;
                PrintAlignmentTable(alignment2);
#endif
                annotate.count_D++;
                annotate.type=2;
            }
            else if( matched_bps_1 >= cutoff_match && alignment1.mismatches <= cutoff_mismatch 
                && alignment1.ref_begin <= cutoff_mismatch && alignment1.query_end >= pBarcode.size()-cutoff_soft_clipping-1 
                && matched_bps_1 > matched_bps_2 )
            {
#ifdef DEBUG
                cout << "P: " << pBarcode << endl;
                PrintAlignmentTable(alignment1);
#endif
                annotate.count_P++;
                annotate.type=1;
            }
            else if( matched_bps_2 >= cutoff_match && alignment2.mismatches <= cutoff_mismatch 
                && alignment2.ref_begin <= cutoff_mismatch && alignment2.query_end >= pBarcode.size()-cutoff_soft_clipping-1 
                && matched_bps_1 < matched_bps_2 )
            {
#ifdef DEBUG
                cout << "D: " << dBarcode << endl;
                PrintAlignmentTable(alignment2);
#endif
                annotate.count_D++;
                annotate.type=2;
            }
            else{
#ifdef DEBUG
                cout <<"No hit" <<endl;
                PrintAlignmentTable(alignment1);
                PrintAlignmentTable(alignment2);
#endif
            }
            return annotate; 
}

/* read exception file */
map<string,int> readExceptionFile(string exception_list){
    map<string,int> map_exception;
    string line;
 
    ifstream infile;
    infile.open(exception_list.c_str());
    if ( ! infile )
        cerr << "fail to open file" << exception_list << endl;

    while( getline( infile, line, '\n' ) ){
        map_exception[line]=1;
    }

    return map_exception;
}

/* check file is exist */
bool is_file_exist(string fileName)
{
    ifstream infile;
    infile.open(fileName.c_str());
    return infile.good();
}
