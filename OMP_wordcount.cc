#include <stdio.h>
#include <queue>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <vector> // for 2D vector
#include <fstream>

using namespace std;
omp_lock_t lock1,lock2, lock3, lock4,lock5; //lock1: MapQ, lock2: files_read, lock3: Htable, lock4: filesQ
    char final_files[20][50] =  {"1.txt","2.txt","3.txt","4.txt","5.txt","6.txt","7.txt","8.txt","9.txt","10.txt","11.txt","12.txt","13.txt","14.txt","15.txt","16.txt","17.txt","18.txt","19.txt","20.txt"};
void init_lock(){
    omp_init_lock(&lock1);
    omp_init_lock(&lock2);
    omp_init_lock(&lock3);
    omp_init_lock(&lock4);
    omp_init_lock(&lock5);
}
void destroy_lock(){
    omp_destroy_lock(&lock1);
    omp_destroy_lock(&lock2);
    omp_destroy_lock(&lock3);
    omp_destroy_lock(&lock4);
    omp_destroy_lock(&lock5);
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
};

int hashFunction(char* x, int buckets) {
    int hashVal = 0, len = strlen(x), i;
    for (i=0; i<len; i++) {
        hashVal += x[i]; 
    }
    return (hashVal % buckets); 
} 

void displayHash(vector<vector<struct hash_item> > &table) { 
    int i,j;
    for (i = 0; i < table.size(); i++) {
        for(j=0; j<table[i].size(); j++)  {
            cout << table[i][j].key<< ":" <<table[i][j].val << "-->";
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
    else 
    {   // given key in hash table. Increment value.
        table[(tid*buckets)+index][i].val += 1;
    }
}

void insertItem2(struct hash_item x, vector<vector<struct hash_item> > &table, int buckets, int idx) {
    int index = idx%buckets, i=0; //idx is index in MHtable, index is index in smaller RHtable.
    for(i=0; i<table[index].size(); i++) {
        if(strcmp(table[index][i].key, x.key)==0) {
            break;
        }
    }
    if (i == table[index].size())    {  //given key not in hash table. Add key.
        struct hash_item new_item;
        strcpy(new_item.key,x.key);
        new_item.val = x.val;
        table[index].push_back(new_item);
    }
    else {                // given key in hash table. Add both values.
        table[index][i].val += x.val;
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

void Reducer(vector<vector<struct hash_item> > &RHtable, vector<vector<struct hash_item> > &MHtable, int buckets, int num_threads) {
    int tid = omp_get_thread_num();
    int i, j, index;
    for (i=0; i<buckets; i++) {
        index = i*num_threads+tid; 
        while(MHtable[index].size() > 0) {
            struct hash_item curr_item;
            curr_item = MHtable[index].back(); MHtable[index].pop_back(); //we don't need lock on MHtable as each reducer will be working on different buckets.
            insertItem2(curr_item, RHtable, buckets, index); //we don't need lock on RHtable as well as each reducer would be adding words to different buckets of RHtable.
        }
    }
}


void displayHash1(vector<vector<struct hash_item> > &table) { 
    int i,j;
    for (i = 0; i < table.size(); i++) {
        for(j=0; j<table[i].size(); j++)  {
            cout << table[i][j].key<< ":" <<table[i][j].val;
        }
        printf("\n");
    }
}
void Writer(vector<vector<struct hash_item> > &RHtable, int buckets, int num_threads) {
    ofstream fout;
    int i, j, index;
    int tid = omp_get_thread_num();
    fout.open(final_files[tid],ios::app); 
    for (i=0; i<buckets/num_threads; i++) {
        index = i*num_threads+tid;
        while(RHtable[index].size() > 0) {
            struct hash_item curr_item;
            curr_item = RHtable[index].back(); RHtable[index].pop_back(); //we don't need lock on MHtable as each reducer will be working on different buckets.
            fout<<curr_item.key<<":"<<curr_item.val<<endl;
        }
    }
    fout.close();
}

int main(int argc, char** argv){
    int maxfiles = atoi(argv[1])*(atoi(argv[3]));
    int num_threads = atoi(argv[2]);
    omp_set_num_threads(num_threads);
    char files[18][50] =  {"1399.txt.utf-8.txt","2600.txt.utf-8.txt","27916.txt.utf-8.txt","3183.txt.utf-8.txt", "36034.txt.utf-8.txt","39290-0.txt","39294.txt.utf-8.txt","39296.txt.utf-8.txt","600.txt.utf-8.txt","2554.txt.utf-8.txt","2638.txt.utf-8.txt", "28054-tei.tei","34114.txt.utf-8.txt","39288.txt.utf-8.txt", "39293-0.txt","39295.txt.utf-8.txt","39297.txt.utf-8.txt","986.txt.utf-8.txt"};
    long MapQ_capacity = 2000; // Size of MapQ
    vector<vector<char> > MapQ; // 2D char vector for keeping mapper work queue.
    MapQ.reserve(MapQ_capacity);
    int buckets=(120/num_threads)*num_threads; //keep #buckets close to 120 for any number of threads.
    vector<vector<struct hash_item> > MHtable(num_threads*buckets); //Mapper HashTable of size #threads*#Uniquebuckets
    vector<vector<struct hash_item> > RHtable(buckets); //Reducer HashTable of size #Uniquebuckts. Reducer does not need to calculate the hash function again as same #unique buckets in both.
    int i = 0,j = 1;
    init_lock();
    queue <int> file_queue; // place index of files  into file queue
    int files_read = 0, reduced = 0,files_written=0;
    int written[num_threads];
    for(i=0;i<num_threads;i++){
        written[i] = 0;
    }
    int read_done=0,map_done=0, reduce_done=0,write_done=0;
    int threshold = num_threads / 2;
    double t1=0.0,t2=0.0;
    for(j=0; j<atoi(argv[3]); j++) { 
    for(i=0; i<atoi(argv[1]); i++)    
    file_queue.push(i); 
    }
    t1 = omp_get_wtime();    

    while((write_done==0)){ 
        #pragma omp parallel
        {
            int curr_thread = omp_get_thread_num();
            if(((curr_thread < threshold) && (files_read < maxfiles)) || ((MapQ.size()==0) && map_done==0 )){ //If MapQ is Empty all threads do Reader work irrespective of threshold.
                 int file_index=-1;                
                 omp_set_lock(&lock4);
                 if(!file_queue.empty()) {file_index = file_queue.front(); file_queue.pop();}
                     omp_unset_lock(&lock4);
                 if(file_index!=-1) {
                     Reader(files[file_index], MapQ);
                     omp_set_lock(&lock2);
                     files_read +=1;
                     omp_unset_lock(&lock2);
                 }
            }
            else if(((curr_thread >= threshold) || (maxfiles==files_read)) && (map_done == 0) ){//if all files are read, all threads should do mapper work irrespective of threshold
                Mapper(MapQ, MHtable, buckets);      
            }
            else if(written[curr_thread]==0) {
/*Implementation1: Too complex: Each Mapper which becomes free (enters this loop) & becomes a Reducer. It would start working on the section of MHtable that it was updating. After that to work on another thread's section of the MHtable, the thread needs to be sure that the other thread has done mapping else it gets wrong count. Also different sections of MHtable have to be marked as read by each thread so that they are not reduced more than once. All this though can be done is too complex to code and would result in a significant overhead in the average case. Thus we would go with implementation 2 which is simpler 
Implementation 2: Simple. Simply put a barrier after the mapper so that each thread becomes reducer simultaneously. No need to keep track of other threads' work. The elseif condition above works as a barrier as threads come out of this loop on map_done=1 and thus explicit barrier instruction need not be put. */
                Reducer(RHtable, MHtable, buckets, num_threads); 
                Writer(RHtable,buckets, num_threads); 
                written[curr_thread] =1;
                omp_set_lock(&lock5);
                files_written+=1;
                omp_unset_lock(&lock5);
            }    
        }
        #pragma omp master
        { //Adjust Threshold
            if ((maxfiles==files_read) && MapQ.empty()) map_done=1;
            if (files_written == num_threads) write_done=1;
            if (MapQ.size() < (0.25* MapQ_capacity) && threshold < (num_threads-1) && (maxfiles != files_read)) //MapQ shrinking rapidly.
              {threshold++; }
            if(MapQ.size() > (0.75 * MapQ_capacity) && threshold > 0) //MapQ increasing rapidly.
              {threshold--; }
        }
    }
    t2 = omp_get_wtime();
    #pragma omp master
    {
        printf("Time on %d threads with %d work is %f\n",num_threads, atoi(argv[3]), t2-t1);
    }
    MapQ.clear();
    RHtable.clear();
    MHtable.clear();
    destroy_lock();
}
