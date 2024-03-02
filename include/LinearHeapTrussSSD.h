#pragma once
#include <fstream>
#include <assert.h>
#include "util.h"
#include "MyFile.h"
#include "edge.h"
using ui = uint32_t;
using ui64 = uint64_t;
using uli = unsigned long long;

#define pb push_back

//linear heap based on the disk
//double link efficiently access data on the disk

class ListLinearHeapTruss {
private:
	ui64 n; // number edges
	
	ui key_cap; // the maximum allowed key value

	ui max_key; // possible max key
	ui min_key; // possible min key

	ui64 *heads; // head of doubly-linked list for a specific weight
    

public:
	std::string m_keys;
	std::string m_pres;
	std::string m_nexts;
	ui64 exist_num;
	ui64 total_io;

	ListLinearHeapTruss(ui64 _n, ui _key_cap, std::string base, std::string supports_dir) {
		n = _n;
		key_cap = _key_cap;
		exist_num = 0;
		total_io = 0;

		heads = nullptr;
		if (heads == nullptr) heads = new ui64[key_cap + 1];
		assert(_key_cap <= key_cap);
		min_key = _key_cap; max_key = 0;
		for (ui i = 0; i <= _key_cap; i++){
			heads[i] = n;
		}
		m_keys = supports_dir;
		m_pres = base+"ListLinear.pres";
		m_nexts = base+"ListLinear.nexts";
	}
	~ListLinearHeapTruss() {
		if (heads != nullptr) {
			delete[] heads;
			heads = nullptr;
		}
	}


	void init(ui64 _n, ui _key_cap) {
		MyReadFile fSup(m_keys);
		fSup.fopen( BUFFERED );
		FILE* fPres = fopen(m_pres.c_str(),"wb");
		FILE* fNexts = fopen(m_nexts.c_str(),"wb");
		ui sup;
		for (ui64 i = 0; i < _n; i++){
			fSup.fseek(i*sizeof(uint32_t));
            fSup.fread(&sup,sizeof(uint32_t));
			insert(i, sup, fPres, fNexts);
		}
		fSup.fclose();
		fclose(fPres);
		fclose(fNexts);
	}
	void init(uint64_t _n, MyReadFile &fOff){
		MyReadFile fSup(m_keys);
		fSup.fopen( BUFFERED );
		FILE* fPres = fopen(m_pres.c_str(),"wb");
		FILE* fNexts = fopen(m_nexts.c_str(),"wb");
		ui sup;
		eid_eid tmp;
		for (ui64 i = 0; i < _n; i++){
			fOff.fseek(i*sizeof(eid_eid));
            fOff.fread(&tmp,sizeof(eid_eid));
			fSup.fseek(tmp.first*sizeof(uint32_t));
            fSup.fread(&sup,sizeof(uint32_t));
			insert(i, sup, fPres, fNexts);
		}
		total_io += _n;
		fSup.fclose();
		fclose(fPres);
		fclose(fNexts);
	}

	// insert (id, key) pair into the data structure
	void insert(ui64 id, ui key, FILE *fPres, FILE *fNexts) {
		if(id >= n)
		{
			printf("id: %lu, n: %lu\n",id,n);
		}
		assert(id < n); 
		assert(key <= key_cap);
		//assert(keys[id] > key_cap);
		exist_num +=1;

		fseek(fPres,id*sizeof(ui64),SEEK_SET);
        fwrite(&n,sizeof(ui64),1,fPres);
		total_io += 1;

		ui64 tmp = heads[key];
		fseek(fNexts,id*sizeof(ui64),SEEK_SET);
        fwrite(&tmp,sizeof(ui64),1,fNexts);
		total_io += 1;

		if (heads[key] != n) {
			fseek(fPres,tmp*sizeof(ui64),SEEK_SET);
        	fwrite(&id,sizeof(ui64),1,fPres);
			total_io += 1;
		}
		heads[key] = id;

		if (key < min_key) min_key = key;
		if (key > max_key) max_key = key;
	}

	void insert(ui64 id, ui key, MyReadFile &fSup, MyReadFile &fPres, MyReadFile &fNexts) {
		assert(id < n); assert(key <= key_cap);
		//assert(keys[id] > key_cap);
		exist_num +=1;
		fSup.fseek(id*sizeof(ui));
        fSup.fwrite(&key,sizeof(ui));

		fPres.fseek(id*sizeof(ui64));
        fPres.fwrite(&n,sizeof(ui64));

		ui64 tmp = heads[key];
		fNexts.fseek(id*sizeof(ui64));
        fNexts.fwrite(&tmp,sizeof(ui64));

		if (heads[key] != n) {
			fPres.fseek(tmp*sizeof(ui64));
        	fPres.fwrite(&id,sizeof(ui64));
		}
		heads[key] = id;

		if (key < min_key) min_key = key;
		if (key > max_key) max_key = key;
	}
	void insert(ui64 id, ui key, MyReadFile &fOff, MyReadFile &fSup, MyReadFile &fPres, MyReadFile &fNexts) {
		assert(id < n); 
		if(key > key_cap){
			printf("key: %u, key_cap: %u\n",key,key_cap);
		}
		assert(key <= key_cap);
		//assert(keys[id] > key_cap);	
		eid_eid tmp_;
		exist_num +=1;

		fOff.fseek(id*sizeof(eid_eid));
        fOff.fread(&tmp_,sizeof(eid_eid));

		fSup.fseek(tmp_.first*sizeof(ui));
        fSup.fwrite(&key,sizeof(ui));
		fSup.fseek(tmp_.second*sizeof(ui));
        fSup.fwrite(&key,sizeof(ui));

		fPres.fseek(id*sizeof(ui64));
        fPres.fwrite(&n,sizeof(ui64));

		ui64 tmp = heads[key];
		fNexts.fseek(id*sizeof(ui64));
        fNexts.fwrite(&tmp,sizeof(ui64));
		if (heads[key] != n) {
			fPres.fseek(tmp*sizeof(ui64));
        	fPres.fwrite(&id,sizeof(ui64));
		}
		heads[key] = id;

		if (key < min_key) min_key = key;
		if (key > max_key) max_key = key;
	}



	// remove a vertex from the data structure
	void remove(ui64 id, ui key, MyReadFile &fPres, MyReadFile &fNexts) {
		ui64 pres_id, nexts_id;
		fPres.fseek(id*sizeof(ui64));
		fPres.fread(&pres_id,sizeof(ui64));

		fNexts.fseek(id*sizeof(ui64));
		fNexts.fread(&nexts_id,sizeof(ui64));
		exist_num--;
		
		if (pres_id == n) {
			heads[key] = nexts_id;
			if (nexts_id != n) 
			{
				fPres.fseek(nexts_id*sizeof(ui64));
				fPres.fwrite(&n,sizeof(ui64));
			}
		}
		else {
			ui64 pid = pres_id;
			fNexts.fseek(pid*sizeof(ui64));
			fNexts.fwrite(&nexts_id,sizeof(ui64));
			if (nexts_id != n) 
			{
				fPres.fseek(nexts_id*sizeof(ui64));
				fPres.fwrite(&pid,sizeof(ui64));
			}
		}
	}

	// pop the (id,key) pair with the minimum key value; return true if success, return false otherwise
	bool pop_min(ui64 &id, ui &key, MyReadFile &fPres, MyReadFile &fNexts) {
		if (empty()) return false;

		id = heads[min_key];
		key = min_key;
		exist_num--;

		ui64 next;
		fNexts.fseek(id*sizeof(ui64));
        fNexts.fread(&next,sizeof(ui64));

		heads[min_key] = next;
		if (heads[min_key] != n) 
		{
			fPres.fseek(heads[min_key]*sizeof(ui64));
        	fPres.fwrite(&n,sizeof(ui64));
		}
		return true;
	}

// 	ui get_n() { return n; }
// 	ui get_key_cap() { return key_cap; }
	ui get_key(ui64 id, MyReadFile &fSup) { 
		ui sup;
		fSup.fseek(id*sizeof(uint32_t));
        fSup.fread(&sup,sizeof(uint32_t));
		return sup; 
	}

	ui get_key(ui64 id, MyReadFile &fSup, MyReadFile &fOff) { 
		ui sup;
		eid_eid tmp;
		fOff.fseek(id*sizeof(eid_eid));
        fOff.fread(&tmp,sizeof(eid_eid));
		
		fSup.fseek(tmp.first*sizeof(uint32_t));
        fSup.fread(&sup,sizeof(uint32_t));
		return sup; 
	}
	bool empty() {
		tighten();
		return min_key > max_key;
	}

// 	ui get_num(ui id){ return count[keys[id]];}
	ui get_minkey() { return min_key; }

// 	void get_ids(std::vector<ui> &ids) {
// 		ids.clear();
// 		tighten();
// 		for (ui i = min_key; i <= max_key; i++) {
// 			for (ui id = heads[i]; id != n; id = nexts[id]) {
// 				ids.pb(id);
// 			}
// 		}
// 	}
// 	void get_ids_of_larger_keys(ui *lst, ui &sz, ui key){
// 		assert(key>=min_key && key <= max_key);

// 		for(ui i = key; i <= max_key; i++){
// 			for(ui id = heads[i]; id!=n; id= nexts[id]){
// 				lst[sz++] = id;
// 			}
// 		}
// 	}
// 	void get_ids_keys(std::vector<ui> &ids, std::vector<ui> &_keys) {
// 		ids.clear(); _keys.clear();
// 		tighten();
// 		for (ui i = min_key; i <= max_key; i++) {
// 			for (ui id = heads[i]; id != n; id = nexts[id]) {
// 				ids.pb(id); _keys.pb(id);
// 			}
// 		}
// 	}



// 	ui size() {
// 		tighten();
// 		ui res = 0;
// 		for (ui i = min_key; i <= max_key; i++) for (ui id = heads[i]; id != n; id = nexts[id]) ++res;
// 		return res;
// 	}

// 	// get the (id,key) pair with the maximum key value; return true if success, return false otherwise
// 	bool get_max(ui &id, ui &key) {
// 		if (empty()) return false;

// 		id = heads[max_key];
// 		key = max_key;
// 		assert(keys[id] == key);
// 		return true;
// 	}

// 	// pop the (id,key) pair with the maximum key value; return true if success, return false otherwise
// 	bool pop_max(ui &id, ui &key) {
// 		if (empty()) return false;

// 		id = heads[max_key];
// 		key = max_key;
// 		assert(keys[id] == key);

// 		heads[max_key] = nexts[id];
// 		if (heads[max_key] != n) pres[heads[max_key]] = n;
// 		return true;
// 	}

// 	// get the (id,key) pair with the minimum key value; return true if success, return false otherwise
// 	bool get_min(ui &id, ui &key) {
// 		if (empty()) return false;

// 		id = heads[min_key];
// 		key = min_key;
// 		assert(keys[id] == key);

// 		return true;
// 	}



// 	// increment the key of vertex id by inc
// 	ui increment(ui id, ui inc = 1) {
// 		assert(keys[id] + inc <= key_cap);

// 		ui new_key = keys[id] + inc;

// 		remove(id);
// 		insert(id, new_key);

// 		return new_key;
// 	}

	// decrement the key of vertex id by dec
	ui decrement(ui sup, ui64 id, MyReadFile &fSup, MyReadFile &fPres, MyReadFile &fNexts, ui dec = 1) {
		ui new_key = sup - dec;

		remove(id,sup,fPres,fNexts);
		insert(id,new_key,fSup,fPres,fNexts);

		return new_key;
	}
	ui decrement(ui sup, ui64 id, MyReadFile &fOff, MyReadFile &fSup, MyReadFile &fPres, MyReadFile &fNexts, ui dec = 1) {
		ui new_key = sup - dec;

		remove(id,sup,fPres,fNexts);
		insert(id,new_key,fOff,fSup,fPres,fNexts);

		return new_key;
	}

private:
	void tighten() {
		while (min_key <= max_key && heads[min_key] == n) ++min_key;
		while (min_key <= max_key && heads[max_key] == n) --max_key;
	}
};