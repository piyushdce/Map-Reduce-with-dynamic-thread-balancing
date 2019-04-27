#include <stdio.h>
#include <queue>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <vector> // for 2D vector
#include <fstream> // for 2D vector

char final_files[16][50] =  {"1.txt","2.txt","3.txt","4.txt","5.txt","6.txt","7.txt","8.txt","9.txt","10.txt","11.txt","12.txt","13.txt","14.txt","15.txt","16.txt"};
using namespace std;
omp_lock_t lock1,lock2, lock3, lock4;
void init_locks(){
    omp_init_lock(&lock1);
    omp_init_lock(&lock2);
    omp_init_lock(&lock3);
    omp_init_lock(&lock4);
}
void dest_locks(){
    omp_destroy_lock(&lock1);
    omp_destroy_lock(&lock2);
    omp_destroy_lock(&lock3);
    omp_destroy_lock(&lock4);
}
void Reader(char *fname_input, vector<vector<char> > &MapQ){
     FILE *input_file;
     char path[] = "./Files/";
     input_file=fopen(strcat(path,fname_input),"r");
     char readchar;
     int count=0;
     int wordCount=0;
     vector<vector<char> > work_item(1);
     while((readchar =fgetc(input_file))!=EOF){
         if(readchar == '|'||readchar == '\n'||readchar == ' ' ||readchar == '?'||readchar == ','||readchar == '.'||readchar == '\"'||
         readchar == '!'||readchar == '@'||readchar == '~'||readchar == '#'||readchar == '$'||readchar == '%'||
         readchar == '%'||readchar == '^'||readchar == '&'||readchar == '*'||readchar == '('||readchar == ')'||
          readchar == '-'||readchar == '_'||readchar == '{'||readchar == '}'||readchar == '['||readchar == ']'|| readchar == ':' ||
         readchar == '<'|| readchar == '+' || readchar == '=' || readchar == '\\' || readchar  == ';' || readchar == '>'||readchar == '\r'||readchar == '/'||readchar == '\''||readchar == '`'
         ){
             if(count>0){
                 wordCount++;
                 if(wordCount == 10000) {
                     omp_set_lock(&lock1);
                     MapQ.push_back(work_item.back());
                     omp_unset_lock(&lock1);
                     work_item[0].clear();
                     wordCount = 0;
                 }
                 else{
                     work_item[0].push_back(' '); //separate words by space so that mapper can identify these words.
                 }
                 count=0;
             }
         }else{
              work_item[0].push_back(readchar);
              count++;
         }
     }
     if(work_item[0].size() !=0) { //push the remaining words to the MapQ.
         omp_set_lock(&lock1);
         MapQ.push_back(work_item.back());
         work_item[0].clear();
         omp_unset_lock(&lock1);
     }
     fclose(input_file); 
}


struct hash_item {
    char key[50];
    int val; 
}hash_item;

int hashFunction(char* x, int buckets) {
    int hashVal = 0, len = strlen(x), i;
    for (i=0; i<len; i++) {
        hashVal += x[i]; 
    }
    return (hashVal % buckets); 
} 
int hashFunctionNode(char* x) {
    int hashVal = 0, len = strlen(x), i;
    for (i=0; i<len; i++) {
        hashVal += x[i]; 
    }
    return hashVal; 
} 

void dispPartialHash(vector<vector<struct hash_item> > &table,int rank,int buckets) { 
     int i,j;
     for (i=rank*buckets; i<(rank+1)*(buckets); i++) {
         for(int j=0; j<table[i].size(); j++)  {
             cout <<table[i][j].key<< ":" <<table[i][j].val<<" ";
         }
     }
}

void insertItem(char *x, vector<vector<struct hash_item> > &table, int buckets) {
    int tid = omp_get_thread_num();
    int index = hashFunction(x,buckets), i=0;
    if (index<0) index=0;
    for(i=0; i<table[(tid*buckets)+index].size(); i++) {
        if(strcmp(table[(tid*buckets)+index][i].key, x)==0) {
            break;
        }
    }    
    if (i == table[(tid*buckets)+index].size()){  //given key not in hash table
        struct hash_item new_item;
        strcpy(new_item.key,x);
        new_item.val = 1;
        table[(tid*buckets)+index].push_back(new_item);
    }   
    else {// given key in hash table. Increment value.
        table[(tid*buckets)+index][i].val += 1;
    }
}
void insertItem2(struct hash_item x,vector<vector<struct hash_item> > &RHtable ,int buckets, int idx,int numP,int rank) {
        int index = idx%buckets, i=0; //idx is index in MHtable, index is index in smaller RHtable.
        int table_chosen =  hashFunctionNode(x.key);
        if (table_chosen<0) table_chosen = 0;
        index = index+(buckets*(table_chosen%numP));
        for(i=0; i<RHtable[index].size(); i++) {
            if(strcmp(RHtable[index][i].key, x.key)==0) {
                break;
            }
        }
        if (i == RHtable[index].size()){  //given key not in hash table. Add key.
            struct hash_item new_item;
            strcpy(new_item.key,x.key);
            new_item.val = x.val;
            RHtable[index].push_back(new_item);
        }
        else{// given key in hash table. Update count
            RHtable[index][i].val += x.val;
        }
}

void Mapper (vector<vector<char> > &MapQ, vector<vector<struct hash_item> > &MHtable, int buckets) {
    vector<char> string_item;
    omp_set_lock(&lock1);
    if(MapQ.size() > 0) {
        string_item = MapQ.back();
        MapQ.pop_back();
    }
    else
    {
        string_item.push_back(' '); 
    }
    omp_unset_lock(&lock1);
    int i, j=0, index;
    char word[50];
    for (i=0; i<string_item.size(); i++) {
        if(string_item[i] != ' ') { 
            word[j] = string_item[i];
            j++;
        }   
        else {
            word[j]='\0';
            j=0;
            insertItem(word, MHtable,buckets);
            strcpy(word,"");
        }
    }
}

void Reducer(vector<vector<struct hash_item> > &RHtable, vector<vector<struct hash_item> > &MHtable, int buckets, int num_threads,int processes,int rank) {
    int tid = omp_get_thread_num();
    int i, j, index;
    char word[50];
    for (i=0; i<buckets; i++) {
        index = i*num_threads+tid; 
        while(MHtable[index].size() > 0) {
            struct hash_item curr_item;
            curr_item = MHtable[index].back(); MHtable[index].pop_back(); 
            insertItem2(curr_item, RHtable,buckets, index,processes,rank);
        }
    }
}


void Writer(vector<vector<struct hash_item> > &RHtable, int buckets, int num_threads,int rank,int numP) {
    ofstream fout;
    fout.open(final_files[rank],ios::app); //every process opens its own file
    int index = 0,i=0;
    int start_idx, end_idx;
    start_idx=rank*buckets; end_idx=(rank+1)*buckets;
    for(i=start_idx; i< end_idx; i=i+numP ) {
        while(RHtable[i].size() > 0) {
            struct hash_item curr_item;
            curr_item = RHtable[i].back(); RHtable[i].pop_back(); 
            fout<<curr_item.key<<":"<<curr_item.val<<endl;
        }
    }
    fout.close();
}


//start of main
int main(int argc, char** argv){
    int rank, p, len, provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided); 
    MPI_Comm_size(MPI_COMM_WORLD, &p);   
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    if(argc !=5){
        if(!rank) printf("Usage: ./exec <#files> <#threads/node> <load factor> <numP>\n");
            return(1); 
    }
    init_locks();
    char files[18][50] =  {"1399.txt.utf-8.txt","2600.txt.utf-8.txt","27916.txt.utf-8.txt","3183.txt.utf-8.txt", "36034.txt.utf-8.txt","39290-0.txt","39294.txt.utf-8.txt","39296.txt.utf-8.txt","600.txt.utf-8.txt","2554.txt.utf-8.txt","2638.txt.utf-8.txt", "28054-tei.tei","34114.txt.utf-8.txt","39288.txt.utf-8.txt", "39293-0.txt","39295.txt.utf-8.txt","39297.txt.utf-8.txt","986.txt.utf-8.txt"}; 
    int maxfiles = atoi(argv[1])*(atoi(argv[3])); 
    int num_threads = atoi(argv[2]); 
    int numP = atoi(argv[4]); 
    omp_set_num_threads(num_threads); 
    long MapQ_capacity = 2000; // Size of MapQ 
    vector<vector<char> > MapQ; // 2D char vector for keeping mapper work queue.//run on each rank (1 node = 1 rank)
    MapQ.reserve(MapQ_capacity); 
    int buckets=(160/num_threads)*num_threads; //keep #buckets close to 120 for any number of threads.
    vector<vector<struct hash_item> > MHtable(num_threads*buckets); //Mapper HashTable of size #threads*#Uniquebuckets
    vector<vector<struct hash_item> > RHtable(buckets*numP); //Mapper HashTable of size #threads*#Uniquebuckets
    int i = 0,j = 1,k=0;
    queue <int> file_queue;// place file inex in queue
    int files_read = 0, reduced = 0;
    int read_done=0,map_done=0, reduce_done=0;
    int threshold = num_threads/2;
    double t1=0.0,t2=0.0;
    for(j=0; j<atoi(argv[3]); j++) { 
        for(i=0; i<atoi(argv[1]); i++)    
            file_queue.push(i); 
    }
    MPI_Datatype STRUCT_DATA_TYPE;
    int blocklens[2];
    MPI_Aint indices[2];
    MPI_Datatype old_types[2];
    /* One value of each type */
    blocklens[0] = 50;
    blocklens[1] = 1; 
    old_types[0] = MPI_CHAR;
    old_types[1] = MPI_INT;
    MPI_Address( &hash_item.key, &indices[0] );
    MPI_Address( &hash_item.val, &indices[1] );
    indices[1] = indices[1] - indices[0];
    indices[0] = 0;
    MPI_Type_struct( 2, blocklens, indices, old_types, &STRUCT_DATA_TYPE );
    MPI_Type_commit( &STRUCT_DATA_TYPE );
    t1 = MPI_Wtime();    
    while((reduce_done==0)){ 
            #pragma omp parallel
            {
                int curr_thread = omp_get_thread_num();
                if(((curr_thread <threshold) && (files_read < maxfiles)) || ((MapQ.size()==0) && map_done==0 )){ //If MapQ is Empty all threads do Reader work irrespective of threshold.
                    int file_index=-1;                
                    omp_set_lock(&lock4);
                    if(!file_queue.empty()) {
                        file_index = file_queue.front(); file_queue.pop();
                    }   
                    omp_unset_lock(&lock4);
                    if(file_index!=-1) {
                        Reader(files[file_index], MapQ);
                        omp_set_lock(&lock2);
                        files_read +=1;
                        omp_unset_lock(&lock2);
                    }
                }
                else if(((curr_thread>=threshold) || (maxfiles==files_read)) && (map_done == 0) ){    //if all files are read, all threads should do mapper work irrespective of threshold
                    Mapper(MapQ, MHtable, buckets);      
                }
                else if(reduce_done==0) {
/*Implementation1: Too complex: Each Mapper which becomes free (enters this loop) & becomes a Reducer. It would start working on the section of MHtable that it was updating. After that to work on another thread's section of the MHtable, the thread needs to be sure that the other thread has done mapping else it gets wrong count. Also different sections of MHtable have to be marked as read by each thread so that they are not reduced more than once. All this though can be done is too complex to code and would result in a significant overhead in the average case. Thus we would go with implementation 2 which is simpler 
Implementation 2: Simple. Simply put a barrier after the mapper so that each thread becomes reducer simultaneously. No need to keep track of other threads' work. The elseif condition above works as a barrier as threads come out of this loop on map_done=1 and thus explicit barrier instruction need not be put. */
                    Reducer(RHtable,MHtable, buckets, num_threads,numP,rank); 
                    omp_set_lock(&lock2);
                    reduced +=1;
                    omp_unset_lock(&lock2);
                }
   
           }
           #pragma omp master
           {
               if ((maxfiles==files_read) && MapQ.empty()) map_done=1;
               if (reduced == num_threads) reduce_done=1;
               if (MapQ.size() < (0.25* MapQ_capacity) && threshold < (num_threads-1) && (maxfiles != files_read)) //MapQ shrinking rapidly.
                   threshold++;
               if(MapQ.size() > (0.75 * MapQ_capacity) && threshold > 0) //MapQ increasing rapidly.
                   threshold--; 
           }
    } // end while


      MPI_Barrier(MPI_COMM_WORLD);   // need barrier because different processes will come out of while at diff time.
      dest_locks();
      t2 = MPI_Wtime();
      // *************** Start MPI communication  **********************/
      int recv_val;
      char recv_key[50];
      MPI_Request reqs[2]; MPI_Status stats[2];
      int sword_cnt, rword_cnt; //send word count and recv word count
      struct hash_item recv_item;
      int p1=-1, p2=-1;
      int start_idx, end_idx, jump, push_idx;

      for (p1=0; p1<numP; p1++) {
          for (p2=0; p2<numP; p2++)  {
              if (p1==p2) {continue;}
                start_idx=(p2*buckets)+p2; end_idx=((p2+1)*buckets); jump=numP; 
                for(j=start_idx; j<end_idx; j+=jump) {
                     sword_cnt = RHtable[j].size();
                     if(rank==p1) {MPI_Isend(&sword_cnt,1,MPI_INT,p2,6,MPI_COMM_WORLD,&reqs[0]); MPI_Waitall(1,reqs,stats);} //send bucket size to p2
                     if(rank==p2) {MPI_Recv(&rword_cnt,1,MPI_INT,p1,6,MPI_COMM_WORLD,&stats[1]);} //recv bucket size from p1
                     if(rank==p1) {
                         MPI_Request request[sword_cnt]; MPI_Status status[sword_cnt];
                         for(int k=0; k<sword_cnt; k++) { //traverse bucket now.
                             MPI_Isend(&RHtable[j][k],1,STRUCT_DATA_TYPE,p2,k,MPI_COMM_WORLD,&request[k]); //send one struct to p2.
                         } // end for k
                         MPI_Waitall(sword_cnt,request,status);
                     } //end if(rank==p1)
                     if(rank==p2) {
                         push_idx = (j)-p2+p1; //push into empty bucket in own region.
                         MPI_Request request[rword_cnt]; MPI_Status status[rword_cnt];
                         for(int k=0; k<rword_cnt; k++) { //recv bucket now.
                             MPI_Recv(&recv_item,1,STRUCT_DATA_TYPE,p1,k,MPI_COMM_WORLD,&status[k]); //recv one struct from p1.
                             RHtable[push_idx].push_back(recv_item); // add the received struct into vector for merging later on.
                         } // end for k
                     } //end if(rank==p2)
                } //end for j
          } //end for p2
     } //end for p1
      // *************** End MPI communication  **********************/

    // ************* Start combining received buckets ***********************/
// Method: Each process does combining for its set of 120 buckets. #process buckets have to be combined into 1 bucket. Final #buckets each process has would be #buckets/#process. Eg if 4 process each process combines 4 buckets into 1 and final #buckets on each process is 30.

    start_idx=rank*buckets; end_idx=(rank+1)*buckets; struct hash_item comb_item;
    for(i=start_idx; i< end_idx; i=i+numP ) {
        for(j=1; j<numP; j++) {
            while(RHtable[(i+j)].size() > 0)    {
                comb_item = RHtable[i+j].back(); RHtable[i+j].pop_back(); 
                insertItem2(comb_item, RHtable, buckets, i, numP,rank);
            }
        }
    }
    // ************* End combining received buckets ***********************/
      double t3 = MPI_Wtime();
      Writer(RHtable,buckets,num_threads,rank,numP);
      if (rank == numP-1) printf("Time(Total:Communication) on rank: %d with process:thd/process:work/process %d:%d:%d is %f:%f \n",rank, numP, num_threads, atoi(argv[3]), t3-t1, t3-t2);
      RHtable.clear();
      MHtable.clear();
      MapQ.clear();
      MPI_Type_free(&STRUCT_DATA_TYPE);
      MPI_Finalize(); 
}
