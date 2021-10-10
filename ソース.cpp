#define _CRT_SECURE_NO_DEPRECATE 1
#include <stdio.h>
#include <math.h>
#include <GL/glut.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <fstream>
#include <random>
using namespace std;

#define N 101
#define WSIZE 600
#define A 3		//運搬車の台数 att48[3],eil101[3],pcb442[5],pr2392[7]
#define aa 1	//コストのa att48[2],eil[1],pcb442[2],pr2392[3]
#define G 800	//遺伝子の数
#define OPT 25	//1.5-opt近傍

void display();
void draw_solution(int rt[N], double position[N][2]);
void keyboard(unsigned char key, int x, int y);
void resize(int w, int h);
void random_route(int rt[N], int seed);
void gravity();
void routeInsert(int temp[], int n, int newdata, int k);
void evABC_sort(double ev_sort[]);
double ev_start2k_k2end(int rt_ABC[], int start_or_end, int  rt);
void render_string(float x, float y, const char* str, double cost);
double cost3(double dist_A, double dist_B, double dist_C);
double cost5(double dist_A, double dist_B, double dist_C, double dist_D, double dist_E);
double dist_ABC(int route[], int count);
void two_opt(int* route, int count);
void insert_opt(int* route, int count);
double insert_but_gr(int route[], int insert_city, int j);
void each_route_min(int route[], int insert_city, double gr_k, int count, double& min, int& min_rt, bool& insert_gr_start, bool& insert_end_rt);
void insert_position(int route[], int& count, bool insert_gr_start, bool insert_end_rt, int insert_city, int min_rt);
void swap(int& x, int& y);
void search_in_each_route(int* rt_1, int* rt_2, int count_1, int count_2);
void two_route_search(int* rt_A, int* rt_B, int* rt_C, int* rt_D, int* rt_E);
//int* rt_F, int* rt_G);
double dist_two_route(int route[], int count, int split_1, int split_2);
double dist_2opt(int rt_temp[], int start, int end);
int get_rand(int min_num,int max_num);

int route[N];		//解（訪問順序）
double pos[N][2];	//町の座標
int num_car = A;

double gr_x;	//デポ（重心）x座標
double gr_y;	//デポ（重心）y座標

int rt_A[N], countA = 0;
int rt_B[N], countB = 0;
int rt_C[N], countC = 0;
int rt_D[N], countD = 0;
int rt_E[N], countE = 0;
//int rt_F[N], countF = 0;
//int rt_G[N], countG = 0;

int best_countA = 0;
int best_countB = 0;
int best_countC = 0;
int best_countD = 0;
int best_countE = 0;
//int best_countF = 0;
//int best_countG = 0;

int best_A[N], best_B[N], best_C[N], best_D[N], best_E[N], best_F[N], best_G[N];
bool rt_zero = true;
double min_cost = 100;
int r_gene[G + 1][N];
int cp_gene[G + 1][N];	//局所最適解から脱出用
int gg = 0;		//遺伝子の数を数えるよう

void idle() {

	random_route(route, rand());

}

int main(int argc, char* argv[])
{
	FILE* fp;
	int i, j;

	srand(time(NULL));
	//グラフィック用関数．削除するな！
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(WSIZE, WSIZE);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutCreateWindow("配送計画問題");
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutKeyboardFunc(keyboard);//s:start; w:wait; q:quit
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glutIdleFunc(idle);

	//ファイルを開く
	if ((fp = fopen("eil101.txt", "r")) == NULL) {
		printf("no file\n");
		exit(0);
	}

	//配達先座標を２次元配列posに読みこむ．
	for (i = 0; i < N; i++) {
		fscanf(fp, "%d,%lf,%lf\n", &j, &(pos[i][0]), &(pos[i][1]));
	}
	fclose(fp);
	glutMainLoop();//グラフィック用関数．削除するな！

}

void random_route(int rt[N], int seed) {
	int v[N], max, cid;
	double gr_A, gr_B, gr_C, A_k, B_k, C_k, gr_k, ev_A, ev_B, ev_C;
	double gr_D, gr_E, gr_F, gr_G, D_k, E_k, F_k, G_k, ev_D, ev_E, ev_F, ev_G, gr_D_end, gr_E_end, gr_F_end, gr_G_end;
	double ev_sort[A] = { 0 };
	double dist_A, dist_B, dist_C, dist_D, dist_E, dist_F, dist_G, now_cost = 0;
	bool insert_end_rtA = false;
	bool insert_end_rtB = false;
	bool insert_end_rtC = false;
	bool insert_end_rtD = false;
	bool insert_end_rtE = false;
	bool insert_end_rtF = false;
	bool insert_end_rtG = false;
	//一次元ずつnewで確保する
	int** taboo_list = new int* [G];
	for (int i = 0; i < G; i++) {
		taboo_list[i] = new int[N];
	}

	bool taboo_one = true;
	int taboo_ct = 0;	//タブーリストの数を数えるため
	int same_ct = 0;	//配列の値が同じ数を数えるため
	bool not_taboo = true;	//タブーリストに含まれているかどうか
	const int same_num = N-1;
	//int kakuritu = 0;
	int ct_gene = 0;	//コピーしたcp_geneのカウントのため
	std::random_device seed_gen;
	std::mt19937 engine(seed_gen());

	if (rt_zero) {
		
		for (int g = 0; g < G; g++) {

			//randam_array(rt);
			for (int i = 0; i < N; i++) {
				std::uint32_t result = rand();
				v[i] = result;
				//cout << v[i] << endl;
			}

			for (int i = 0; i < N; i++) {
				max = -1;
				for (int j = 0; j < N; j++) {
					if (v[j] > max) {
						max = v[j];
						cid = j;
					}
				}
				rt[i] = cid;
				v[cid] = -1;

			}


			//最初だけ実行される
			if (taboo_one) {
				//cout << "一度だけ" << endl;
				taboo_one = false;

				for (int p = 0; p < same_num; p++) {
					taboo_list[taboo_ct][p] = rt[p];
					//cout << "taboo_list[" << p << "]:" << taboo_list[g][p] << endl;
				}
				taboo_ct++;
				for (int k = 0; k < N; k++) {
					r_gene[g][k] = rt[k];
					//cout << "r_gene[" << k << "]:" << r_gene[g][k] << endl << endl;
				}
			}
			else {
				//r_geneに追加
				/*for (int k = 0; k < N; k++) {
					r_gene[g][k] = rt[k];
					//cout << "r_gene[" << g << "][" << k << "]:" << r_gene[g][k] << endl << endl;
				}*/

				for (int i = 0; i < taboo_ct; i++) {

					for (int j = 0; j < same_num; j++) {
						for (int k = 0; k < same_num; k++) {

							//タブーリストに値が含まれているかどうか
							if (rt[j] == taboo_list[i][k]) {
								same_ct++;
								//cout << "same_ct:" << same_ct << endl;
							}
						}
					}
					if (same_ct == same_num) {
						//if (same_start<=same_ct&&same_ct<=same_end) {
						not_taboo = false;
						//cout << "同じ" << endl;
					}
					same_ct = 0;
				}

				//タブーリストに含まれている間は探索する

				if (!not_taboo) {
					//cout << "実行" << endl;
					for (int i = 0; i < N; i++) {
						std::uint32_t result = rand();
						v[i] = result;
						//cout << v[i] << endl;
					}

					for (int i = 0; i < N; i++) {
						max = -1;
						for (int j = 0; j < N; j++) {
							if (v[j] > max) {
								max = v[j];
								cid = j;
							}
						}
						rt[i] = cid;
						v[cid] = -1;
					}
					for (int i = 0; i < taboo_ct; i++) {

						for (int j = 0; j < same_num; j++) {
							for (int k = 0; k < same_num; k++) {

								//タブーリストに値が含まれているかどうか
								if (rt[j] == taboo_list[i][k]) {
									same_ct++;
									//cout << "same_ct:" << same_ct << endl;
								}
							}
						}
						if (same_ct != same_num) {
							not_taboo = true;
							//cout << "同じ" << endl;
						}
						same_ct = 0;
					}
					not_taboo = true;
				}

				//タブーリストに含まれていないものがr_geneとなる
				for (int k = 0; k < N; k++) {
					r_gene[g][k] = rt[k];
					//cp_gene[g][k] = rt[k];
					//cout << "r_gene[" << g << "][" << k << "]:" << r_gene[g][k] << endl << endl;
				}
				
				//タブーリストに追加

				for (int p = 0; p < same_num; p++) {
					taboo_list[taboo_ct][p] = rt[p];

				}
				taboo_ct++;

				not_taboo = true;
			}
		}
		cout << taboo_ct << endl;
		rt_zero = false;
	}

	//一次元ずつ解放する
	for (int i = 0; i < G; i++) {
		delete[] taboo_list[i];
	}
	delete[] taboo_list;

	countA = 0;
	countB = 0;
	countC = 0;
	countD = 0;
	countE = 0;
	int mutation_ct = 0;;
	int swap1 = 0, swap2 = 0;
	int element = 0;	//コピーした配列の要素数を数えるため

	if (gg < G) {

		rt_A[countA] = r_gene[gg][0];
		rt_B[countB] = r_gene[gg][1];
		rt_C[countC] = r_gene[gg][2];
		/*rt_D[countD] = r_gene[gg][3];
		rt_E[countE] = r_gene[gg][4];
		*/

		gravity();
		//cout << gr_x << ":" << gr_y << endl;

		gr_A = sqrt(pow((gr_x - pos[rt_A[0]][0]), 2) + pow((gr_y - pos[rt_A[0]][1]), 2));
		gr_B = sqrt(pow((gr_x - pos[rt_B[0]][0]), 2) + pow((gr_y - pos[rt_B[0]][1]), 2));
		gr_C = sqrt(pow((gr_x - pos[rt_C[0]][0]), 2) + pow((gr_y - pos[rt_C[0]][1]), 2));
		/*gr_D = sqrt(pow((gr_x - pos[rt_D[0]][0]), 2) + pow((gr_y - pos[rt_D[0]][1]), 2));
		gr_E = sqrt(pow((gr_x - pos[rt_E[0]][0]), 2) + pow((gr_y - pos[rt_E[0]][1]), 2));
		*/

		A_k = sqrt(pow((pos[rt_A[0]][0] - pos[r_gene[gg][A]][0]), 2) + pow((pos[rt_A[0]][1] - pos[r_gene[gg][A]][1]), 2));
		B_k = sqrt(pow((pos[rt_B[0]][0] - pos[r_gene[gg][A]][0]), 2) + pow((pos[rt_B[0]][1] - pos[r_gene[gg][A]][1]), 2));
		C_k = sqrt(pow((pos[rt_C[0]][0] - pos[r_gene[gg][A]][0]), 2) + pow((pos[rt_C[0]][1] - pos[r_gene[gg][A]][1]), 2));
		/*D_k = sqrt(pow((pos[rt_D[0]][0] - pos[r_gene[gg][A]][0]), 2) + pow((pos[rt_D[0]][1] - pos[r_gene[gg][A]][1]), 2));
		E_k = sqrt(pow((pos[rt_E[0]][0] - pos[r_gene[gg][A]][0]), 2) + pow((pos[rt_E[0]][1] - pos[r_gene[gg][A]][1]), 2));
		*/

		gr_k = sqrt(pow((gr_x - pos[r_gene[gg][A]][0]), 2) + pow((gr_y - pos[r_gene[gg][A]][1]), 2));


		ev_A = A_k + gr_k - gr_A;
		ev_B = B_k + gr_k - gr_B;
		ev_C = C_k + gr_k - gr_C;
		/*ev_D = D_k + gr_k - gr_D;
		ev_E = E_k + gr_k - gr_E;
		*/

		ev_sort[0] = ev_A;
		ev_sort[1] = ev_B;
		ev_sort[2] = ev_C;
		/*ev_sort[3] = ev_D;
		ev_sort[4] = ev_E;
		*/

		evABC_sort(ev_sort);

		//挿入法の結果、一番距離が短いものがどれか
		//一番距離が短いものルートに都市を追加する
		if (ev_A == ev_sort[0]) {
			countA++;
			rt_A[countA] = r_gene[gg][A];

		}
		else if (ev_B == ev_sort[0]) {
			countB++;
			rt_B[countB] = r_gene[gg][A];

		}
		else if (ev_C == ev_sort[0]) {
			countC++;
			rt_C[countC] = r_gene[gg][A];

		}
		/*else if (ev_D == ev_sort[0]) {
			countD++;
			rt_D[countD] = r_gene[gg][A];
		}
		else if (ev_E == ev_sort[0]) {
			countE++;
			rt_E[countE] = r_gene[gg][A];
		}
		*/

		//残りの都市がどのルートに入ればいいか調べる
		//今は、ルートにそれぞれ１つは都市が入っている＆１つのルートだけ２つ目の都市が入っている
		for (int i = A + 1; i < N; i++) {
			gr_k = sqrt(pow((gr_x - pos[r_gene[gg][i]][0]), 2) + pow((gr_y - pos[r_gene[gg][i]][1]), 2));

			//ルートAだけを考える
			double min_A = ev_start2k_k2end(rt_A, 0, r_gene[gg][i]) + gr_k;
			int min_rt_A = rt_A[0];//デポの間（デポと最初の都市の間）に挿入
			bool insert_gr_startA = true;
			each_route_min(rt_A, i, gr_k, countA, min_A, min_rt_A, insert_gr_startA, insert_end_rtA);


			//ルートBだけ考える
			double min_B = ev_start2k_k2end(rt_B, 0, r_gene[gg][i]) + gr_k;
			int min_rt_B = rt_B[0];//デポの間（デポと最初の都市の間）に挿入
			bool insert_gr_startB = true;
			each_route_min(rt_B, i, gr_k, countB, min_B, min_rt_B, insert_gr_startB, insert_end_rtB);

			//ルートCだけ考える
			double min_C = ev_start2k_k2end(rt_C, 0, r_gene[gg][i]) + gr_k;
			int min_rt_C = rt_C[0];//デポの間（デポと最初の都市の間）に挿入
			bool insert_gr_startC = true;
			each_route_min(rt_C, i, gr_k, countC, min_C, min_rt_C, insert_gr_startC, insert_end_rtC);

			//ルートDだけ考える
			/*double min_D = ev_start2k_k2end(rt_D, 0, r_gene[gg][i]) + gr_k;
			int min_rt_D = rt_D[0];//デポの間（デポと最初の都市の間）に挿入
			bool insert_gr_startD = true;
			each_route_min(rt_D, i, gr_k, countD, min_D, min_rt_D, insert_gr_startD, insert_end_rtD);
			//ルートEだけ考える
			double min_E = ev_start2k_k2end(rt_E, 0, r_gene[gg][i]) + gr_k;
			int min_rt_E = rt_E[0];//デポの間（デポと最初の都市の間）に挿入
			bool insert_gr_startE = true;
			each_route_min(rt_E, i, gr_k, countE, min_E, min_rt_E, insert_gr_startE, insert_end_rtE);
			*/

			//A・B・Cの評価を並べる
			ev_sort[0] = min_A;
			ev_sort[1] = min_B;
			ev_sort[2] = min_C;
			/*ev_sort[3] = min_D;
			ev_sort[4] = min_E;
			*/

			evABC_sort(ev_sort);


			if (min_A == ev_sort[0]) {

				insert_position(rt_A, countA, insert_gr_startA, insert_end_rtA, r_gene[gg][i], min_rt_A);

			}
			else if (min_B == ev_sort[0]) {
				insert_position(rt_B, countB, insert_gr_startB, insert_end_rtB, r_gene[gg][i], min_rt_B);

			}
			else if (min_C == ev_sort[0]) {
				insert_position(rt_C, countC, insert_gr_startC, insert_end_rtC, r_gene[gg][i], min_rt_C);

			}
			/*else if (min_D == ev_sort[0]) {
				insert_position(rt_D, countD, insert_gr_startD, insert_end_rtD, r_gene[gg][i], min_rt_D);
			}
			else if (min_E == ev_sort[0]) {
				insert_position(rt_E, countE, insert_gr_startE, insert_end_rtE, r_gene[gg][i], min_rt_E);
			}
			*/

			insert_end_rtA = false;
			insert_end_rtB = false;
			insert_end_rtC = false;
			insert_end_rtD = false;
			insert_end_rtE = false;
			
			insert_gr_startA = false;
			insert_gr_startB = false;
			insert_gr_startC = false;
			/*insert_gr_startD = false;
			insert_gr_startE = false;
			*/
		}
		
		//局所最適解からの脱出のために初期解をコピーする
		for (int c = 0; c < countA+1; c++) {
			cp_gene[gg][c] = rt_A[c];
			ct_gene++;
		}
		element = ct_gene;
		for (int c = 0; c < countB+1; c++) {
			cp_gene[gg][c + element] = rt_B[c];
			ct_gene++;
		}
		element = ct_gene;
		for (int c = 0; c < countC+1; c++) {
			cp_gene[gg][c +element] = rt_C[c];
			ct_gene++;
		}
		element = ct_gene;
		/*for (int c = 0; c < countD+1; c++) {
			cp_gene[gg][c + element] = rt_D[c];
			ct_gene++;
		}
		element = ct_gene;
		for (int c = 0; c < countE+1; c++) {
			cp_gene[gg][c + element] = rt_E[c];
			ct_gene++;
		}
		*/
		
		//1つの経路ごとの1.5opt近傍
		if (countA > 1) {
			//cout << "A突然変異" << endl;
			insert_opt(rt_A, countA);//突然変異
		}

		if (countB > 1) {
			//cout << "B突然変異" << endl;
			insert_opt(rt_B, countB);//突然変異
		}

		if (countC > 1) {
			//cout << "C突然変異" << endl;
			insert_opt(rt_C, countC);//突然変異
		}

		/*if (countD > 1) {
			//cout << "C突然変異" << endl;
			insert_opt(rt_D, countD);//突然変異
		}
		if (countE > 1) {
			//cout << "C突然変異" << endl;
			insert_opt(rt_E, countE);//突然変異
		}
		*/
		
		//1つの経路ごとの2opt近傍
		
		if (countA > 4) {
			//cout << "A突然変異" << endl;
			two_opt(rt_A, countA);//突然変異
		}

		if (countB > 4) {
			//cout << "B突然変異" << endl;
			two_opt(rt_B, countB);//突然変異
		}

		if (countC > 4) {
			//cout << "C突然変異" << endl;
			two_opt(rt_C, countC);//突然変異
		}

		/*if (countD > 4) {
			//cout << "C突然変異" << endl;
			two_opt(rt_D, countD);//突然変異
		}
		if (countE > 4) {
			//cout << "C突然変異" << endl;
			two_opt(rt_E, countE);//突然変異
		}
		*/
	

		//2つの経路ごとの2opt近傍
		
		if (countA > 1 && countB > 1 && countC > 1) {
			//&& countD > 1 && countE && countF > 1 && countG > 1) {
			two_route_search(rt_A, rt_B, rt_C, rt_D, rt_E);
		}
				

		//今のコスト
		dist_A = dist_ABC(rt_A, countA);
		dist_B = dist_ABC(rt_B, countB);
		dist_C = dist_ABC(rt_C, countC);
		dist_D = dist_ABC(rt_D, countD);
		dist_E = dist_ABC(rt_E, countE);

		if (A == 3) {
			now_cost = cost3(dist_A, dist_B, dist_C);
		}
		else if (A == 5) {
			now_cost = cost5(dist_A, dist_B, dist_C, dist_D, dist_E);
		}

		if (now_cost < min_cost) {
			cout << "---------------------更新" << endl;
			min_cost = now_cost;

			for (int p = 0; p < countA + 1; p++) {
				best_A[p] = rt_A[p];
			}
			for (int p = 0; p < countB + 1; p++) {
				best_B[p] = rt_B[p];
			}
			for (int p = 0; p < countC + 1; p++) {
				best_C[p] = rt_C[p];
			}
			/*for (int p = 0; p < countD + 1; p++) {
				best_D[p] = rt_D[p];
			}
			for (int p = 0; p < countE + 1; p++) {
				best_E[p] = rt_E[p];
			}
			*/

			best_countA = countA;
			best_countB = countB;
			best_countC = countC;
			/*best_countD = countD;
			best_countE = countE;
			*/

		}
		
		//局所最適解の脱出-------------------------------------------------------------ここから
		
		//step1 突然変異
		mutation_ct = get_rand(1, 25);	//突然変異が実行される確率
		//cout << "突然変異実行回数:" << mutation_ct << endl;
		for (int m = 0; m < mutation_ct; m++) {
			swap1 = get_rand(0,N-1);
			swap2 = get_rand(0,N-1);
			//cout << "swap1:" << swap1 << "  swap2:" << swap2 << endl;
			swap(cp_gene[gg][swap1], cp_gene[gg][swap2]);
		}
		ct_gene = 0;

		//step2　コピー初期解を作成
		for (int p = 0; p < countA + 1; p++) {
			rt_A[p] = cp_gene[gg][p];
			ct_gene++;
		}
		element = ct_gene;
		for (int p = 0; p < countB + 1; p++) {
			rt_B[p] = cp_gene[gg][p+element];
			ct_gene++;
		}
		element = ct_gene;
		for (int p = 0; p < countC + 1; p++) {
			rt_C[p] = cp_gene[gg][p + element];
			ct_gene++;
		}
		element = ct_gene;
		/*for (int p = 0; p < countD + 1; p++) {
			rt_D[p] = cp_gene[gg][p + element];
			ct_gene++;
		}
		element = ct_gene;
		for (int p = 0; p < countE + 1; p++) {
			rt_E[p] = cp_gene[gg][p + element];
			ct_gene++;
		}
		*/
		element = 0;

		//step3
		//1つの経路ごとの1.5opt近傍
		if (countA > 1) {
			//cout << "A突然変異" << endl;
			insert_opt(rt_A, countA);//突然変異
		}

		if (countB > 1) {
			//cout << "B突然変異" << endl;
			insert_opt(rt_B, countB);//突然変異
		}

		if (countC > 1) {
			//cout << "C突然変異" << endl;
			insert_opt(rt_C, countC);//突然変異
		}

		/*if (countD > 1) {
			//cout << "C突然変異" << endl;
			insert_opt(rt_D, countD);//突然変異
		}
		if (countE > 1) {
			//cout << "C突然変異" << endl;
			insert_opt(rt_E, countE);//突然変異
		}
		*/

		//1つの経路ごとの2opt近傍

		if (countA > 4) {
			//cout << "A突然変異" << endl;
			two_opt(rt_A, countA);//突然変異
		}

		if (countB > 4) {
			//cout << "B突然変異" << endl;
			two_opt(rt_B, countB);//突然変異
		}

		if (countC > 4) {
			//cout << "C突然変異" << endl;
			two_opt(rt_C, countC);//突然変異
		}

		/*if (countD > 4) {
			//cout << "C突然変異" << endl;
			two_opt(rt_D, countD);//突然変異
		}
		if (countE > 4) {
			//cout << "C突然変異" << endl;
			two_opt(rt_E, countE);//突然変異
		}
		*/

		//2つの経路ごとの2opt近傍

		/*if (countA > 1 && countB > 1 && countC > 1) {
			//&& countD > 1 && countE) {
			two_route_search(rt_A, rt_B, rt_C, rt_D, rt_E);
		}*/

		//step4　コストを求めて更新するか決める
		//今のコスト
		dist_A = dist_ABC(rt_A, countA);
		dist_B = dist_ABC(rt_B, countB);
		dist_C = dist_ABC(rt_C, countC);
		dist_D = dist_ABC(rt_D, countD);
		dist_E = dist_ABC(rt_E, countE);

		if (A == 3) {
			now_cost = cost3(dist_A, dist_B, dist_C);
		}
		else if (A == 5) {
			now_cost = cost5(dist_A, dist_B, dist_C, dist_D, dist_E);
		}

		if (now_cost < min_cost) {
			cout << "---------------------脱出更新" << endl;
			min_cost = now_cost;

			for (int p = 0; p < countA + 1; p++) {
				best_A[p] = rt_A[p];
			}
			for (int p = 0; p < countB + 1; p++) {
				best_B[p] = rt_B[p];
			}
			for (int p = 0; p < countC + 1; p++) {
				best_C[p] = rt_C[p];
			}
			/*for (int p = 0; p < countD + 1; p++) {
				best_D[p] = rt_D[p];
			}
			for (int p = 0; p < countE + 1; p++) {
				best_E[p] = rt_E[p];
			}
			*/

			best_countA = countA;
			best_countB = countB;
			best_countC = countC;
			/*best_countD = countD;
			best_countE = countE;
			*/

		}

		//-----------------------------------------------------------------ここまで

		cout << gg << endl;
		gg++;
		if (gg == G) {
			ofstream writing_file;
			string filename = "result.csv";
			writing_file.open(filename, std::ios::app);
			double writing_num = min_cost;
			writing_file << writing_num << endl;
			writing_file.close();
		}

	}


	glutPostRedisplay();

}



void draw_solution(int rt[N], double position[N][2]) {

	int i;
	char string[100];

	glClear(GL_COLOR_BUFFER_BIT);

	glColor3d(0.2, 0.2, 0.2);
	glBegin(GL_LINE_LOOP);
	glVertex2d(0.00, 0.00);
	glVertex2d(0.00, 1.00);
	glVertex2d(1.00, 1.00);
	glVertex2d(1.00, 0.00);
	glEnd();

	//-----rt_A-----
	glColor3d(1.0, 1.0, 0.0); //yellow
	glBegin(GL_LINE_LOOP);
	//glVertex2d(gr_x, gr_y)

	for (int i = 0; i < best_countA + 1; i++) {
		glVertex2dv(position[best_A[i]]);
	}
	glVertex2d(gr_x, gr_y);
	glEnd();

	//-----rt_B-----
	glColor3d(1.0, 0.0, 0.0); //red
	glBegin(GL_LINE_LOOP);
	//glVertex2d(gr_x, gr_y);
	for (int i = 0; i < best_countB + 1; i++) {
		glVertex2dv(position[best_B[i]]);
	}
	glVertex2d(gr_x, gr_y);
	glEnd();

	//-----rt_C-----
	glColor3d(1.0, 1.0, 1.0); //white
	glBegin(GL_LINE_LOOP);
	//glVertex2d(gr_x, gr_y);
	for (int i = 0; i < best_countC + 1; i++) {
		glVertex2dv(position[best_C[i]]);
	}
	glVertex2d(gr_x, gr_y);
	glEnd();

	//-----rt_D-----
	/*glColor3d(1.0, 0.0, 1.0); //white
	glBegin(GL_LINE_LOOP);
	//glVertex2d(gr_x, gr_y);
	for (int i = 0; i < best_countD + 1; i++) {
		glVertex2dv(position[best_D[i]]);
	}
	glVertex2d(gr_x, gr_y);
	glEnd();

	//-----rt_E-----
	glColor3d(0.0, 0.0, 1.0); //white
	glBegin(GL_LINE_LOOP);
	//glVertex2d(gr_x, gr_y);
	for (int i = 0; i < best_countE + 1; i++) {
		glVertex2dv(position[best_E[i]]);
	}
	glVertex2d(gr_x, gr_y);
	glEnd();
	//-----rt_F-----
	glColor3d(0.0, 1.0, 1.0); //white
	glBegin(GL_LINE_LOOP);
	//glVertex2d(gr_x, gr_y);
	for (int i = 0; i < best_countF + 1; i++) {
		glVertex2dv(position[best_F[i]]);
	}
	glVertex2d(gr_x, gr_y);
	glEnd();
	//-----rt_G-----
	glColor3d(0.0, 1.0, 0.0); //white
	glBegin(GL_LINE_LOOP);
	//glVertex2d(gr_x, gr_y);
	for (int i = 0; i < best_countG + 1; i++) {
		glVertex2dv(position[best_G[i]]);
	}
	glVertex2d(gr_x, gr_y);
	glEnd();*/

	glColor3d(0.0, 0.0, 1.0); //blue
	glPointSize(3);
	glBegin(GL_POINTS);
	for (i = 0; i < N; i++) {
		glVertex2dv(position[i]);
	}
	glEnd();

	//デポの表示
	glColor3d(1.0, 0.0, 0.0);//red
	glPointSize(5);
	glBegin(GL_POINTS);

	glVertex2d(gr_x, gr_y);
	glEnd();

	render_string(0.68f, 1.02f, "cost=", min_cost);

	glutSwapBuffers();
}


void display() {
	draw_solution(route, pos);
}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 's':
		glutIdleFunc(idle);
		break;
	case 'w':
		glutPostRedisplay();
		glutIdleFunc(0);
		break;
	case 'q':
		exit(1);
		break;
	}
}

void resize(int w, int h) {
	double margin = 0.1;

	glViewport(0, 0, w, h);
	glLoadIdentity();
	glOrtho(0.0 - margin, 1.0 + margin, 0.0 - margin, 1.0 + margin, -1.0, 1.0);
	glOrtho(-w / (WSIZE * 1.0), w / (WSIZE * 1.0), -h / (WSIZE * 1.0), h / (WSIZE * 1.0), -1.0, 1.0);
}


void drawBitmapString(void* font, char* string)
{

	glPushAttrib(GL_CURRENT_BIT);

	/* ビットマップ文字列の描画 */
	while (*string)
		glutBitmapCharacter(font, *string++);

	glPopAttrib();
}

void gravity() {
	double x_all = 0, y_all = 0;

	for (int i = 0; i < N; i++) {
		x_all += pos[i][0];
		y_all += pos[i][1];
	}

	gr_x = x_all / N;
	gr_y = y_all / N;

}

int get_rand(int min_num,int max_num) {

	std::random_device rnd;     // 非決定的な乱数生成器を生成
	std::mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
	std::uniform_int_distribution<> rand100(min_num, max_num);        


	return rand100(mt);
}

void routeInsert(int temp[], int n, int newdata, int k) {

	for (int i = n; i > k; i--) {
		temp[i] = temp[i - 1];
	}
	temp[k] = newdata;
}

void evABC_sort(double ev_sort[]) {

	for (int s = 0; s < A - 1; s++) {
		for (int t = s + 1; t < A; t++) {
			if (ev_sort[t] < ev_sort[s]) {
				double temp = ev_sort[t];
				ev_sort[t] = ev_sort[s];
				ev_sort[s] = temp;
			}
		}
	}
}

double ev_start2k_k2end(int rt_ABC[], int start_or_end, int rt) {
	double gr_ABC, ABC_k;

	gr_ABC = sqrt(pow((gr_x - pos[rt_ABC[start_or_end]][0]), 2) + pow((gr_y - pos[rt_ABC[start_or_end]][1]), 2));

	ABC_k = sqrt(pow((pos[rt_ABC[start_or_end]][0] - pos[rt][0]), 2) + pow((pos[rt_ABC[start_or_end]][1] - pos[rt][1]), 2));

	return ABC_k - gr_ABC;

}

void render_string(float x, float y, const char* str, double cost) { //文字表記用
	float z = -1.0f;
	glColor3f(1, 1, 1);
	glRasterPos3d(x, y, z);
	char* c = (char*)str;

	std::string num = std::to_string(cost);
	const char* costnum = num.c_str();

	while (*c != '\0') {
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *c++);
	}

	while (*costnum != '\0') {
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *costnum++);
	}
}

double cost3(double dist_A, double dist_B, double dist_C) {
	double sum = 0, ave = 0;

	sum = dist_A + dist_B + dist_C;
	ave = sum / A;

	return sum + aa * (fabs(ave - dist_A) + fabs(ave - dist_B) + fabs(ave - dist_C)) / A;

}

double cost5(double dist_A, double dist_B, double dist_C, double dist_D, double dist_E) {
	double sum = 0, ave = 0;

	sum = dist_A + dist_B + dist_C + dist_D + dist_E;
	ave = sum / A;

	return sum + aa * (fabs(ave - dist_A) + fabs(ave - dist_B) + fabs(ave - dist_C) + fabs(ave - dist_D) + fabs(ave - dist_E)) / A;

}

//１つのルートの距離
double dist_ABC(int route[], int count) {
	double distance = 0;

	distance += sqrt(pow((gr_x - pos[route[0]][0]), 2) + pow((gr_y - pos[route[0]][1]), 2));
	for (int i = 0; i < count; i++) {
		distance += sqrt(pow((pos[route[i]][0] - pos[route[i + 1]][0]), 2) + pow((pos[route[i]][1] - pos[route[i + 1]][1]), 2));
	}
	distance += sqrt(pow((gr_x - pos[route[count]][0]), 2) + pow((gr_y - pos[route[count]][1]), 2));

	return distance;
}

//2つのルートどうしでの距離
double dist_two_route(int route[], int num, int split_1, int split_2) {
	double distance = 0;

	//デポとの距離を求める
	//distance += sqrt(pow((gr_x - pos[route[split_1-1]][0]), 2) + pow((gr_y - pos[route[split_1-1]][1]), 2));
	//distance += sqrt(pow((gr_x - pos[route[split_2-1]][0]), 2) + pow((gr_y - pos[route[split_2 - 1]][1]), 2));
	for (int i = 0; i < num - 1; i++) {
		distance += sqrt(pow((pos[route[i]][0] - pos[route[i + 1]][0]), 2) + pow((pos[route[i]][1] - pos[route[i + 1]][1]), 2));
	}

	return distance;
}

double dist_2opt(int rt_temp[], int start, int end) {
	double distance = 0;
	//逆順する最初の都市の1つ前の都市との距離
	distance += sqrt(pow((pos[route[start - 1]][0] - pos[route[start]][0]), 2) + pow((pos[route[start - 1]][1] - pos[route[start]][1]), 2));
	//逆順する最後の都市の1つ後の都市との距離
	distance += sqrt(pow((pos[route[end - 1]][0] - pos[route[end]][0]), 2) + pow((pos[route[end - 1]][1] - pos[route[end]][1]), 2));

	for (int i = 0; i < end - start; i++) {
		distance += sqrt(pow((pos[rt_temp[i]][0] - pos[rt_temp[i + 1]][0]), 2) + pow((pos[rt_temp[i]][1] - pos[rt_temp[i + 1]][1]), 2));
	}

	return distance;
}

void insert_opt(int* route, int count) {
	//1.5-opt近傍
	bool is_improved = true;//改善されたか
	double dist_ac, dist_cb, dist_ab, now_dist;
	double min_mutation = 100;//変化量が最小なものを保持
	double now_mutation = 0;//今の変化量
	int insert_city = 0;//挿入する都市
	int pre_city = 0;//この都市の次に挿入
	double min_dist = dist_ABC(route, count);//局所探索する前の距離

	while (is_improved) {
		for (int k = 0; k < count; k++) {
			for (int c = 0; c < k; c++) {
				dist_ab = sqrt(pow((pos[route[k]][0] - pos[route[k + 1]][0]), 2) + pow((pos[route[k]][1] - pos[route[k + 1]][1]), 2));
				dist_ac = sqrt(pow((pos[route[k]][0] - pos[route[c]][0]), 2) + pow((pos[route[k]][1] - pos[route[c]][1]), 2));
				dist_cb = sqrt(pow((pos[route[c]][0] - pos[route[k + 1]][0]), 2) + pow((pos[route[c]][1] - pos[route[k + 1]][1]), 2));

				now_mutation = dist_ac + dist_cb - dist_ab;

				if (now_mutation < min_mutation) {
					min_mutation = now_mutation;
					insert_city = c;
					pre_city = k;
				}
			}
			for (int c = k + 2; c < count + 1; c++) {
				dist_ab = sqrt(pow((pos[route[k]][0] - pos[route[k + 1]][0]), 2) + pow((pos[route[k]][1] - pos[route[k + 1]][1]), 2));
				dist_ac = sqrt(pow((pos[route[k]][0] - pos[route[c]][0]), 2) + pow((pos[route[k]][1] - pos[route[c]][1]), 2));
				dist_cb = sqrt(pow((pos[route[c]][0] - pos[route[k + 1]][0]), 2) + pow((pos[route[c]][1] - pos[route[k + 1]][1]), 2));

				now_mutation = dist_ac + dist_cb - dist_ab;

				if (now_mutation < min_mutation) {
					min_mutation = now_mutation;
					insert_city = c;
					pre_city = k;
				}
			}

		}

		if (insert_city > pre_city) {

			for (int s = insert_city; s > pre_city + 1; s--) {
				swap(route[s], route[s - 1]);//前の配列と交換
			}
		}
		else {
			for (int s = insert_city; s < pre_city; s++) {
				swap(route[s], route[s + 1]);//次の配列と交換
			}
		}

		now_dist = dist_ABC(route, count);

		if (now_dist < min_dist) {
			min_dist = now_dist;
			//cout << "1.5-opt近傍実行" << endl;
			is_improved = true;
		}
		else {
			//もとに戻す
			if (insert_city > pre_city) {
				for (int s = pre_city + 1; s < insert_city; s++) {
					swap(route[s], route[s + 1]);//前の配列と交換
				}
			}
			else {
				for (int s = insert_city; s < pre_city; s++) {
					swap(route[s], route[s + 1]);//前の配列と交換
				}
			}
			is_improved = false;
		}
	}

}

void two_opt(int* route, int count) {
	double dist_ac, dist_cb, dist_ab, now_dist;// , dist_pn, dist_ab, dist_pc, dist_cn, now_dist;
	//double min_mutation = 100;//変化量が最小なものを保持
	//double now_mutation = 0;//今の変化量
	//int insert_city = 0;//挿入する都市
	//int pre_city = 0;//この都市の次に挿入
	//bool is_improved = true;//改善されたか
	bool is_improved_2opt = true;
	double min_dist = dist_ABC(route, count);//局所探索する前の距離

	//2-opt近傍
	int* now_temp = new int[N];
	int min_x = 0, min_c = 0;
	//値をコピー
	for (int x = 0; x < count + 1; x++) {
		now_temp[x] = route[x];
	}
	
	double pre_dist = dist_ABC(route, count);//局所探索法（2-opt近傍）の前の距離

	while (is_improved_2opt) {
		for (int x = 0; x < count; x++) {
			for (int c = x + 1; c < count + 1; c++) {

				for (int k = 0; k <= c - x; k++) {
					now_temp[x + k] = route[c - k];
				}

				now_dist = dist_ABC(now_temp, count);

				if (now_dist < min_dist) {
					min_dist = now_dist;
					min_x = x;
					min_c = c;
				}
				//元に戻す
				for (int x = 0; x < count + 1; x++) {
					now_temp[x] = route[x];
				}
			}
		}

		if (min_dist < pre_dist) {
			pre_dist = min_dist;
			//経路を更新
			//cout << "2-opt実行" << endl;
			for (int k = 0; k <= min_c - min_x; k++) {
				now_temp[min_x + k] = route[min_c - k];
			}
			for (int k = min_x; k <= min_c; k++) {
				route[k] = now_temp[k];
			}

		}
		else {
			is_improved_2opt = false;
		}

	}

	//1.5-opt近傍
	/*while (is_improved) {
		for (int k = 0; k < count; k++) {
			for (int c = 0; c < k; c++) {
				dist_ab = sqrt(pow((pos[route[k]][0] - pos[route[k + 1]][0]), 2) + pow((pos[route[k]][1] - pos[route[k + 1]][1]), 2));
				dist_ac = sqrt(pow((pos[route[k]][0] - pos[route[c]][0]), 2) + pow((pos[route[k]][1] - pos[route[c]][1]), 2));
				dist_cb = sqrt(pow((pos[route[c]][0] - pos[route[k + 1]][0]), 2) + pow((pos[route[c]][1] - pos[route[k + 1]][1]), 2));

				now_mutation = dist_ac + dist_cb - dist_ab;

				if (now_mutation < min_mutation) {
					min_mutation = now_mutation;
					insert_city = c;
					pre_city = k;
				}
			}
			for (int c = k + 2; c < count + 1; c++) {
				dist_ab = sqrt(pow((pos[route[k]][0] - pos[route[k + 1]][0]), 2) + pow((pos[route[k]][1] - pos[route[k + 1]][1]), 2));
				dist_ac = sqrt(pow((pos[route[k]][0] - pos[route[c]][0]), 2) + pow((pos[route[k]][1] - pos[route[c]][1]), 2));
				dist_cb = sqrt(pow((pos[route[c]][0] - pos[route[k + 1]][0]), 2) + pow((pos[route[c]][1] - pos[route[k + 1]][1]), 2));

				now_mutation = dist_ac + dist_cb - dist_ab;

				if (now_mutation < min_mutation) {
					min_mutation = now_mutation;
					insert_city = c;
					pre_city = k;
				}
			}

		}

		if (insert_city > pre_city) {

			for (int s = insert_city; s > pre_city + 1; s--) {
				swap(route[s], route[s - 1]);//前の配列と交換
			}
		}
		else {
			for (int s = insert_city; s < pre_city; s++) {
				swap(route[s], route[s + 1]);//次の配列と交換
			}
		}

		now_dist = dist_ABC(route, count);

		if (now_dist < min_dist) {
			min_dist = now_dist;
			//cout << "1.5-opt近傍実行" << endl;
			is_improved = true;
		}
		else {
			//もとに戻す
			if (insert_city > pre_city) {
				for (int s = pre_city + 1; s < insert_city; s++) {
					swap(route[s], route[s + 1]);//前の配列と交換
				}
			}
			else {
				for (int s = insert_city; s < pre_city; s++) {
					swap(route[s], route[s + 1]);//前の配列と交換
				}
			}
			is_improved = false;
		}
	}*/

	


	delete[] now_temp;

}



void swap(int& x, int& y) {
	int temp;    // 値を一時保存する変数

	temp = x;
	x = y;
	y = temp;
}



double insert_but_gr(int route[], int insert_city, int j) {
	double X_Y = 0, X_k = 0, Y_k = 0;

	X_Y = sqrt(pow((pos[route[j]][0] - pos[route[j + 1]][0]), 2) + pow((pos[route[j]][1] - pos[route[j + 1]][1]), 2));
	X_k = sqrt(pow((pos[route[j]][0] - pos[r_gene[gg][insert_city]][0]), 2) + pow((pos[route[j]][1] - pos[r_gene[gg][insert_city]][1]), 2));
	Y_k = sqrt(pow((pos[route[j + 1]][0] - pos[r_gene[gg][insert_city]][0]), 2) + pow((pos[route[j + 1]][1] - pos[r_gene[gg][insert_city]][1]), 2));

	return X_k + Y_k - X_Y;
}

void each_route_min(int route[], int insert_city, double gr_k, int count, double& min, int& min_rt, bool& insert_gr_start, bool& insert_end_rt) {

	double ev_XY, ev_end;


	if (count != 0) {
		for (int j = 0; j < count; j++) {
			//ルート内で挿入（デポ以外）
			ev_XY = insert_but_gr(route, insert_city, j);

			if (min > ev_XY) {
				min = ev_XY;
				min_rt = j;//この都市の次にｋを挿入
				insert_gr_start = false;

			}
		}

		//デポに帰ってくる前に挿入
		ev_end = ev_start2k_k2end(route, count, r_gene[gg][insert_city]) + gr_k;

		if (min > ev_end) {
			min = ev_end;

			insert_end_rt = true;
			insert_gr_start = false;

		}

	}


}

void insert_position(int route[], int& count, bool insert_gr_start, bool insert_end_rt, int insert_city, int min_rt) {
	//デポと最初の都市の間に挿入
	if (insert_gr_start) {
		count++;
		routeInsert(route, count, insert_city, 0);

	}
	//最後とデポの間に挿入
	else if (insert_end_rt) {
		count++;
		route[count] = insert_city;

	}
	else {
		count++;
		routeInsert(route, count, insert_city, min_rt + 1);

	}
}

//2-optを行う2つのルートの組合せ
void two_route_search(int* rt_A, int* rt_B, int* rt_C, int* rt_D, int* rt_E) {

	if (num_car == 3) {
		//att48[3]のときに実行

		search_in_each_route(rt_A, rt_B, countA, countB);

		search_in_each_route(rt_B, rt_C, countB, countC);
		search_in_each_route(rt_C, rt_A, countC, countA);

	}
	else if (num_car == 5) {
		//pcb442[5]のときに実行
		search_in_each_route(rt_A, rt_B, countA, countB);
		search_in_each_route(rt_B, rt_C, countB, countC);
		search_in_each_route(rt_C, rt_D, countC, countD);
		search_in_each_route(rt_D, rt_E, countD, countE);
		search_in_each_route(rt_E, rt_A, countE, countA);
	}
	
}

void search_in_each_route(int* rt_1, int* rt_2, int count_1, int count_2) {
	int split_1 = count_1 + 1;		//ルートを分けるとき用
	int split_2 = count_2 + 1;	//ルートを分けるとき用
	int num = split_1 + split_2;	//繋げたときの配列の数
	//int newtemp[N] = { 0 };
	int* newtemp = new int[num];
	int start_2opt_position = 0;	//2optをする範囲のスタート
	int end_2opt_position = 0;		//2optをする範囲のエンド

	double min_dist = 0, now_dist = 0;
	double rt_1_to_rt_2 = 0, rt_2_to_rt_1 = 0;

	//rt_1の最後とrt_2の最初をつなげるべきか
	rt_1_to_rt_2 = sqrt(pow((pos[rt_1[count_1]][0] - pos[rt_2[0]][0]), 2) + pow((pos[rt_1[count_1]][1] - pos[rt_2[0]][1]), 2));

	//rt_2の最後とrt_1の最初をつなげるべきか
	rt_2_to_rt_1 = sqrt(pow((pos[rt_2[count_2]][0] - pos[rt_1[0]][0]), 2) + pow((pos[rt_2[count_2]][1] - pos[rt_1[0]][1]), 2));

	//3つのルートを一つのルートへ
	if (rt_1_to_rt_2 < rt_2_to_rt_1) {
		for (int i = 0; i < split_1; i++) {
			newtemp[i] = rt_1[i];
		}
		for (int i = split_1; i < num; i++) {
			newtemp[i] = rt_2[i - split_1];
		}
	}
	else {
		for (int i = 0; i < split_2; i++) {
			newtemp[i] = rt_2[i];
		}
		for (int i = split_2; i < num; i++) {
			newtemp[i] = rt_1[i - split_2];
		}
	}


	//現在の距離（2つのルートの距離）を最短とする
	min_dist = dist_two_route(newtemp, num, split_1, split_2);


	//2opt近傍の実行
	for (int s = 0; s <= count_1; s++) {
		for (int e = count_2; e > count_1; e--) {

			for (int i = s, j = e; i < count_1 + 1; i++, j--) {
				swap(newtemp[i], newtemp[j]);
			}

			//実行後の距離を求める
			now_dist = dist_two_route(newtemp, num, split_1, split_2);

			//最短にであるか
			if (now_dist < min_dist) {
				//cout << "22opt" << endl;
				//最短であるときの2optをする範囲のはじめとおわりを保持する
				start_2opt_position = s;
				end_2opt_position = e;
			}

			//逆順を元に戻す
			for (int i = s, j = e; i < count_1 + 1; i++, j--) {
				swap(newtemp[i], newtemp[j]);
			}
			now_dist = 0;
		}

	}

	//元のルートを最短のルートに変更
	//最短が更新されたときに2opt実行
	if (start_2opt_position != 0 && end_2opt_position != 0) {
		//cout << "実行・・・・・・・・・" << endl;
		for (int i = start_2opt_position, j = end_2opt_position; i < count_1 + 1; i++, j--) {
			swap(newtemp[i], newtemp[j]);
		}

		//まとめたルートを2つのルートに戻し、最短に更新
		if (rt_1_to_rt_2 < rt_2_to_rt_1) {
			for (int i = 0; i < split_1; i++) {
				rt_1[i] = newtemp[i];
			}
			for (int j = split_1; j < num; j++) {
				rt_2[j - split_1] = newtemp[j];
			}
		}
		else {
			for (int i = 0; i < split_2; i++) {
				rt_2[i] = newtemp[i];
			}
			for (int x = 0; x < split_1; x++) {
				rt_1[x] = newtemp[x + split_2];
			}
		}

	}
	delete[] newtemp;
}

