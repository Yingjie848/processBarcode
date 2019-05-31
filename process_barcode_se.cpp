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
    string seq_exception_list; 
    string output_prefix="out";

    int c;
    while ((c=getopt(argc,argv,"m:M:s:a:p:")) != -1)
    {
        switch (c)
        {
            case 'm' : if(optarg) cutoff_match = atoi(optarg); break; 
            case 'M' : if(optarg) cutoff_mismatch = atoi(optarg); break; 
            case 's' : if(optarg) cutoff_soft_clipping = atoi(optarg); break; 
            case 'a' : if(optarg) seq_exception_list=optarg; break;
            case 'p' : if(optarg) output_prefix=optarg; break;
        }
    }
    string seq_file = argv[optind++];

    string out_seq_mapped_reads = output_prefix+".map.fa";
    string out_seq_xhoI_reads = output_prefix+".map.xhoI.fa";
    string out_seq_telo_reads = output_prefix+".map.telo.fa";
    string out_SPB = output_prefix+".unmap.SPB.fa"; //Output Sequence labeled with Proximal Barcode (SPB)
    string out_SDB = output_prefix+".unmap.SDB.fa"; //Output Sequence labeled with Distal Barcode (SDB)
    string out_stat = output_prefix+".stat";

    //output stat file
    ofstream stat ( out_stat.c_str() );
    if ( ! stat )
        cerr << "fail to create stat file" << out_stat <<endl;

    stat << "#Parameters: " << " -m " << cutoff_match << " -M " << cutoff_mismatch 
         << " -a " << seq_exception_list << " -p " << output_prefix 
         << " R1.fa: " << seq_file ; 
    
    // Declare
    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment1;
    StripedSmithWaterman::Alignment alignment2;
    // Declares map
    map<string,int> map_seq_exception;
    map<string,int> map_pBarcode; 
    map<string,int>::iterator map_pBarcode_Iter; 
    map<string,int> map_dBarcode; 
    map<string,int>::iterator map_dBarcode_Iter; 
    map<string,int> map_pattern; 
    map<string,int>::iterator map_pattern_Iter; 
    //
    string seq;  //input sequence
    int count_seq=0;       //count total # of reads
    int count_seq_exception=0; //count # of excepted reads
    int count_seq_xhoI=0;      //count xhoI
    int count_seq_telo=0;      //count telomere
    int count_seq_unmapped=0; //count unmapped reads
    int count_seq_P=0;    //count R1 proximal reads
    int count_seq_D=0;    //count R1 distal reads
    int count_seq_usable_proximal_barcode=0;  //count usable proximal barcode, include P-P, P-D, D-P,P-None, None-P
    int count_seq_usable_distal_barcode=0;  //count usable distal barcode, include P-D, D-P

    //open files
    ifstream infile;
    infile.open(seq_file.c_str());
    if ( ! infile )
        cerr << "fail to open file" << seq_file << endl;
     
    //output files
    ofstream map_seq ( out_seq_mapped_reads.c_str() );
    if ( ! map_seq )
        cerr << "fail to create file" << out_seq_mapped_reads <<endl;

    ofstream map_seq_xhoI ( out_seq_xhoI_reads.c_str() );
    if ( ! map_seq_xhoI )
        cerr << "fail to create file" << out_seq_xhoI_reads <<endl;

    ofstream map_seq_telo ( out_seq_telo_reads.c_str() );
    if ( ! map_seq_telo )
        cerr << "fail to create file" << out_seq_telo_reads <<endl;

    ofstream clean1 ( out_SPB.c_str() );
    if ( ! clean1 )
        cerr << "fail to create file" << out_SPB <<endl;
    
    ofstream clean2 ( out_SDB.c_str() );
    if ( ! clean2 )
        cerr << "fail to create file" << out_SDB <<endl;

    //load exception list
    if( is_file_exist(seq_exception_list) )
        map_seq_exception=readExceptionFile(seq_exception_list);


    while( getline( infile, seq, '\n' ) ) //read line by '\n'
    {
            count_seq++; //count total sequence
            if( map_seq_exception[boost::lexical_cast<string>(count_seq-1)]==1 )
            {
                //cout << count_seq << " " <<  map_seq_exception[boost::lexical_cast<string>(count_seq-1)] << endl;
                count_seq_exception++;
                map_seq << seq << endl;

                //align xhoI
                if( seq.find(xhoI) == 0 ){
                    count_seq_xhoI++;
                    map_seq_xhoI << seq << endl;
                }else if( (seq.find(telomere1) == 0) || (seq.find(telomere2) == 0) || (seq.find(telomere3) == 0) || (seq.find(telomere4) == 0) || (seq.find(telomere5) == 0) || (seq.find(telomere6) == 0) ){
                    count_seq_telo++;
                    map_seq_telo << seq << endl;
                } 
                continue;
            }
#ifdef DEBUG
            cout << "#" << count_seq <<endl; 
#endif            
            count_seq_unmapped++;
            /*     process seq     */
            //align proximal barcode, search barcode using string::find function, then do alignment using smith-waterman method
            alignment1 = align_query_to_reference(pBarcode,sub_pBarcode,seq);
            int matched_bps_1=(alignment1.ref_end-alignment1.ref_begin+1)-alignment1.mismatches;
            float matched_perc_1=(float)matched_bps_1*100/(alignment1.ref_end-alignment1.ref_begin+1);
            
            //align distal barcode
            alignment2 = align_query_to_reference(dBarcode,sub_dBarcode,seq); 
            int matched_bps_2=(alignment2.ref_end-alignment2.ref_begin+1)-alignment2.mismatches;
            float matched_perc_2=(float)matched_bps_2*100/(alignment2.ref_end-alignment2.ref_begin+1);
#ifdef DEBUG
            //cout << seq << endl;
            cout << "P:" << matched_bps_1 << "/" << matched_perc_1 << "%\tD:" << matched_bps_2 << "/" << matched_perc_2 << "%" << endl;
#endif

            //annotate
            //filter by matched base pairs, mismatch number, reference start of alignment, query end of alignment, percent of matched base pairs
            ANNOTATE_T annotate = annotate_seq(cutoff_match, cutoff_mismatch, cutoff_soft_clipping, matched_bps_1, matched_bps_2, matched_perc_1, matched_perc_2, alignment1, alignment2);
            count_seq_P+=annotate.count_P;
            count_seq_D+=annotate.count_D;

            //output clean sequence
            if( annotate.type == 1 )      //clean seq P
            {
                map_pBarcode[seq.substr(alignment1.ref_begin,alignment1.ref_end-alignment1.ref_begin+1)]++;    
                clean1 << seq.substr(alignment1.ref_end+1,seq.size()-alignment1.ref_end+1) << endl;
            }
            else if( annotate.type == 2 ) //clean seq D
            {
                map_dBarcode[seq.substr(alignment2.ref_begin,alignment2.ref_end-alignment2.ref_begin+1)]++;    
                clean2 << seq.substr(alignment2.ref_end+1,seq.size()-alignment2.ref_end+1) << endl;
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
    //count annotated sequence
    stat << "Total_reads\tException\tXhoI_reads\tTelomere_reads\tTotal_reads_unmapped\tProximal_barcode_annotated\tDistal_barcode_annotated\n"; 
    stat << count_seq << "\t" << count_seq_exception << "\t" << count_seq_xhoI << "\t" << count_seq_telo << "\t" << count_seq_unmapped << "\t" << count_seq_P << "\t" << count_seq_D << endl;

    infile.close();
    clean1.close();clean2.close();

    return 0;
}

/* ==========USAGE function=========== */
static void usage ()
{	
    cout << "Usage: process_barcode -m 9 -M 1 -s 1 -a <exception_list> -p <output_prefix>  <read.fa>\n"
         << " -m <int>      minimum of match number, default is 9\n"
         << " -M <int>      max mismatch number, default is 1\n"
         << " -s <int>      max soft clipping base pairs, default is 1\n"
         << " -a <string>,  Exception list generated by mapping raw reads to genome, only unique ID name in this list\n"
         << " -p <string>,  output prefix for all output files\n"
         << " <read.fa>     input sequence line by line\n";
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
                && alignment2.ref_begin <= cutoff_mismatch && alignment2.query_end >= dBarcode.size()-cutoff_soft_clipping-1 
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
                && alignment2.ref_begin <= cutoff_mismatch && alignment2.query_end >= dBarcode.size()-cutoff_soft_clipping-1 
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
