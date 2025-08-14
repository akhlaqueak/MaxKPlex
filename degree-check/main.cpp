#include "Utility.h"

int main(int argc, char *argv[]) {
		FILE *f = Utility::open_file(dir.c_str(), "rb");
	
		ui tt, n, m, max_degree=0;
		fread(&tt, sizeof(int), 1, f);
		fread(&n, sizeof(int), 1, f);
		fread(&m, sizeof(int), 1, f);
	
		
		ui *degree = new ui[n];
		ui* pstart = new ept[n+1];
		ui* edges = new ui[m];
		
		fread(degree, sizeof(int), n, f);
		pstart[0] = 0;
		for(ui i = 0;i < n;i ++) {
			if(degree[i] > 0) {
				fread(edges+pstart[i], sizeof(int), degree[i], f);
				
				// remove self loops and parallel edges
				ui *buff = edges+pstart[i];
				sort(buff, buff+degree[i]);
				ui idx = 0;
				for(ui j = 0;j < degree[i];j ++) {
					if(buff[j] >= n) printf("vertex id %u wrong\n", buff[j]);
					if(buff[j] == i||(j > 0&&buff[j] == buff[j-1])) continue;
					buff[idx ++] = buff[j];
				}
				degree[i] = idx;
			}
			
			pstart[i+1] = pstart[i] + degree[i];
			max_degree = max(degree[i], max_degree);
		}
		printf("\tn = %s, m = %s (undirected), max_degree = %s \n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str(), Utility::integer_to_string(max_degree).c_str());
		
		fclose(f);
		delete[] degree;
		delete[] pstart;
		delete[] edges;
	}