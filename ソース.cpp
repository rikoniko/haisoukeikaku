#define _CRT_SECURE_NO_DEPRECATE 1
#include <stdio.h>
#include <stdlib.h> //RAND_MAX
#include <math.h>
#include <GL/glut.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
using namespace std;

#define N 101
#define WSIZE 600
#define A 3		//�^���Ԃ̑䐔 att48[3],eil101[3],pcb442[5],pr2392[7]
#define aa 1	//�R�X�g��a att48[2],eil[1],pcb442[2],pr2392[3]
#define G 800	//��`�q�̐�
#define OPT 25	//1.5-opt�ߖT

void display();
void draw_solution(int rt[N], double position[N][2]);
void keyboard(unsigned char key, int x, int y);
void resize(int w, int h);
void random_route(int rt[N], int seed);
void gravity();
void routeInsert(int temp[], int n, int newdata, int k);
void evABC_sort(double ev_sort[]);
double ev_start2k_k2end(int rt_ABC[],int start_or_end, int  rt);
void render_string(float x, float y, const char* str, double cost);
double cost3(double dist_A, double dist_B, double dist_C);
double cost5(double dist_A, double dist_B, double dist_C, double dist_D, double dist_E);
double cost7(double dist_A, double dist_B, double dist_C, double dist_D, double dist_E, double dist_F, double dist_G);
double dist_ABC(int route[], int count);
void mutation(int *route, int count);
int GetRandom(int min, int max);
double insert_but_gr(int route[],int insert_city,int j);
void each_route_min(int route[], int insert_city, double gr_k, int count, double& min, int& min_rt, bool& insert_gr_start, bool& insert_end_rt);
void insert_position(int route[], int &count, bool insert_gr_start, bool insert_end_rt, int insert_city, int min_rt);
void swap(int& x, int& y);
void search_in_each_route(int* rt_1, int* rt_2, int count_1, int count_2);
void two_route_search(int* rt_A, int* rt_B, int* rt_C, int* rt_D, int* rt_E, int* rt_F, int* rt_G);
double dist_two_route(int route[], int count);

int route[N];		//���i�K�⏇���j
double pos[N][2];	//���̍��W
int num_car = A;

double gr_x;	//�f�|�i�d�S�jx���W
double gr_y;	//�f�|�i�d�S�jy���W

int rt_A[N], countA = 0;
int rt_B[N], countB = 0;
int rt_C[N], countC = 0;
int rt_D[N], countD = 0;
int rt_E[N], countE = 0;
int rt_F[N], countF = 0;
int rt_G[N], countG = 0;

int best_countA = 0;
int best_countB = 0;
int best_countC = 0;
int best_countD = 0;
int best_countE = 0;
int best_countF = 0;
int best_countG = 0;

int best_A[N], best_B[N], best_C[N], best_D[N], best_E[N], best_F[N], best_G[N];
bool rt_zero = true;
double min_cost=100;
int r_gene[G + 1][N];
int gg = 0;		//��`�q�̐��𐔂���悤

void idle() {
	random_route(route, rand());
}

void main(int argc, char* argv[])
{
	FILE* fp;
	int i, j;

	srand(time(NULL));

	//�O���t�B�b�N�p�֐��D�폜����ȁI
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(WSIZE, WSIZE);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutCreateWindow("TSP");
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutKeyboardFunc(keyboard);//s:start; w:wait; q:quit
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glutIdleFunc(idle);

	//�t�@�C�����J��
	if ((fp = fopen("eil101.txt", "r")) == NULL) {
		printf("no file\n");
		exit(0);
	}

	//�z�B����W���Q�����z��pos�ɓǂ݂��ށD
	for (i = 0; i < N; i++) {
		fscanf(fp, "%d,%lf,%lf\n", &j, &(pos[i][0]), &(pos[i][1]));
	}
	fclose(fp);

	glutMainLoop();//�O���t�B�b�N�p�֐��D�폜����ȁI
}

void random_route(int rt[N], int seed) {
	int i, j, v[N], max, cid;
	double gr_A, gr_B, gr_C, A_k, B_k, C_k, gr_k, ev_A, ev_B, ev_C;
	double gr_D, gr_E, gr_F,gr_G, D_k, E_k, F_k, G_k, ev_D, ev_E, ev_F,ev_G, gr_D_end, gr_E_end, gr_F_end,gr_G_end;
	double ev_sort[A] = { 0 };
	double dist_A,dist_B,dist_C, dist_D, dist_E, dist_F, dist_G,now_cost=0;
	bool insert_end_rtA = false;
	bool insert_end_rtB = false;
	bool insert_end_rtC = false;
	bool insert_end_rtD = false;
	bool insert_end_rtE = false;
	bool insert_end_rtF = false;
	bool insert_end_rtG = false;

	

	if (rt_zero) {
		for (int g = 0; g < G; g++) {
			for (i = 0; i < N; i++) {
				v[i] = rand();
			}

			for (i = 0; i < N; i++) {
				max = -1;
				for (j = 0; j < N; j++) {
					if (v[j] > max) {
						max = v[j];
						cid = j;
					}
				}
				rt[i] = cid;
				v[cid] = -1;
			}
			for (int k = 0; k < N; k++) {
				r_gene[g][k] = route[k];
			}
		}

		rt_zero = false;
	}
	countA = 0;
	countB = 0;
	countC = 0;
	countD = 0;
	countE = 0;
	countF = 0;
	countG = 0;

	if (gg < G) {
		
		rt_A[countA] = r_gene[gg][0];
		rt_B[countB] = r_gene[gg][1];
		rt_C[countC] = r_gene[gg][2];
		/*rt_D[countD] = r_gene[gg][3];
		rt_E[countE] = r_gene[gg][4];
		rt_F[countF] = r_gene[gg][5];
		rt_G[countG] = r_gene[gg][6];*/

		gravity();
		//cout << gr_x << ":" << gr_y << endl;
		
		gr_A = sqrt(pow((gr_x - pos[rt_A[0]][0]), 2) + pow((gr_y - pos[rt_A[0]][1]), 2));
		gr_B = sqrt(pow((gr_x - pos[rt_B[0]][0]), 2) + pow((gr_y - pos[rt_B[0]][1]), 2));
		gr_C = sqrt(pow((gr_x - pos[rt_C[0]][0]), 2) + pow((gr_y - pos[rt_C[0]][1]), 2));
		/*gr_D = sqrt(pow((gr_x - pos[rt_D[0]][0]), 2) + pow((gr_y - pos[rt_D[0]][1]), 2));
		gr_E = sqrt(pow((gr_x - pos[rt_E[0]][0]), 2) + pow((gr_y - pos[rt_E[0]][1]), 2));
		gr_F = sqrt(pow((gr_x - pos[rt_F[0]][0]), 2) + pow((gr_y - pos[rt_F[0]][1]), 2));
		gr_G = sqrt(pow((gr_x - pos[rt_G[0]][0]), 2) + pow((gr_y - pos[rt_G[0]][1]), 2));*/
		
	
		A_k = sqrt(pow((pos[rt_A[0]][0] - pos[r_gene[gg][A]][0]), 2) + pow((pos[rt_A[0]][1] - pos[r_gene[gg][A]][1]), 2));
		B_k = sqrt(pow((pos[rt_B[0]][0] - pos[r_gene[gg][A]][0]), 2) + pow((pos[rt_B[0]][1] - pos[r_gene[gg][A]][1]), 2));
		C_k = sqrt(pow((pos[rt_C[0]][0] - pos[r_gene[gg][A]][0]), 2) + pow((pos[rt_C[0]][1] - pos[r_gene[gg][A]][1]), 2));
		/*D_k = sqrt(pow((pos[rt_D[0]][0] - pos[r_gene[gg][A]][0]), 2) + pow((pos[rt_D[0]][1] - pos[r_gene[gg][A]][1]), 2));
		E_k = sqrt(pow((pos[rt_E[0]][0] - pos[r_gene[gg][A]][0]), 2) + pow((pos[rt_E[0]][1] - pos[r_gene[gg][A]][1]), 2));
		F_k = sqrt(pow((pos[rt_F[0]][0] - pos[r_gene[gg][A]][0]), 2) + pow((pos[rt_F[0]][1] - pos[r_gene[gg][A]][1]), 2));
		G_k = sqrt(pow((pos[rt_G[0]][0] - pos[r_gene[gg][A]][0]), 2) + pow((pos[rt_G[0]][1] - pos[r_gene[gg][A]][1]), 2));*/

		gr_k = sqrt(pow((gr_x - pos[r_gene[gg][A]][0]), 2) + pow((gr_y - pos[r_gene[gg][A]][1]), 2));


		ev_A = A_k + gr_k - gr_A;
		ev_B = B_k + gr_k - gr_B;
		ev_C = C_k + gr_k - gr_C;
		/*ev_D = D_k + gr_k - gr_D;
		ev_E = E_k + gr_k - gr_E;
		ev_F = F_k + gr_k - gr_F;
		ev_G = G_k + gr_k - gr_G;*/

		ev_sort[0] = ev_A;
		ev_sort[1] = ev_B;
		ev_sort[2] = ev_C;
		/*ev_sort[3] = ev_D;
		ev_sort[4] = ev_E;
		ev_sort[5] = ev_F;
		ev_sort[6] = ev_G;*/

		evABC_sort(ev_sort);

		//�^���ԂR�ƂT�ňႤ��
		if (ev_A == ev_sort[0]) {
			countA++;
			rt_A[countA] = r_gene[gg][A];

		}
		else if (ev_B == ev_sort[0]) {
			countB++;
			rt_B[countB] = r_gene[gg][A];

		}
		else if(ev_C==ev_sort[0]){
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
		else if (ev_F == ev_sort[0]) {
			countF++;
			rt_F[countF] = r_gene[gg][A];

		}
		else if (ev_G == ev_sort[0]) {
			countG++;
			rt_G[countG] = r_gene[gg][A];

		}*/


		for (int i = A+1; i < N; i++) {
			gr_k = sqrt(pow((gr_x - pos[r_gene[gg][i]][0]), 2) + pow((gr_y - pos[r_gene[gg][i]][1]), 2));

			//���[�gA�������l����
			double min_A = ev_start2k_k2end(rt_A, 0,r_gene[gg][i]) + gr_k;
			int min_rt_A = rt_A[0];//�f�|�̊ԁi�f�|�ƍŏ��̓s�s�̊ԁj�ɑ}��
			bool insert_gr_startA = true;
			each_route_min(rt_A, i, gr_k, countA, min_A, min_rt_A, insert_gr_startA, insert_end_rtA);


			//���[�gB�����l����
			double min_B = ev_start2k_k2end(rt_B,0, r_gene[gg][i]) + gr_k;
			int min_rt_B = rt_B[0];//�f�|�̊ԁi�f�|�ƍŏ��̓s�s�̊ԁj�ɑ}��
			bool insert_gr_startB = true;
			each_route_min(rt_B, i, gr_k, countB, min_B, min_rt_B, insert_gr_startB, insert_end_rtB);

			//���[�gC�����l����
			double min_C = ev_start2k_k2end(rt_C,0, r_gene[gg][i]) + gr_k;
			int min_rt_C = rt_C[0];//�f�|�̊ԁi�f�|�ƍŏ��̓s�s�̊ԁj�ɑ}��
			bool insert_gr_startC = true;
			each_route_min(rt_C,i,gr_k,countC,min_C, min_rt_C,insert_gr_startC,insert_end_rtC);		
			
			//���[�gD�����l����
			/*double min_D = ev_start2k_k2end(rt_D, 0, r_gene[gg][i]) + gr_k;
			int min_rt_D = rt_D[0];//�f�|�̊ԁi�f�|�ƍŏ��̓s�s�̊ԁj�ɑ}��
			bool insert_gr_startD = true;
			each_route_min(rt_D, i, gr_k, countD, min_D, min_rt_D, insert_gr_startD, insert_end_rtD);

			//���[�gE�����l����
			double min_E = ev_start2k_k2end(rt_E, 0, r_gene[gg][i]) + gr_k;
			int min_rt_E = rt_E[0];//�f�|�̊ԁi�f�|�ƍŏ��̓s�s�̊ԁj�ɑ}��
			bool insert_gr_startE = true;
			each_route_min(rt_E, i, gr_k, countE, min_E, min_rt_E, insert_gr_startE, insert_end_rtE);

			//���[�gF�����l����
			double min_F = ev_start2k_k2end(rt_F, 0, r_gene[gg][i]) + gr_k;
			int min_rt_F = rt_F[0];//�f�|�̊ԁi�f�|�ƍŏ��̓s�s�̊ԁj�ɑ}��
			bool insert_gr_startF = true;
			each_route_min(rt_F, i, gr_k, countF, min_F, min_rt_F, insert_gr_startF, insert_end_rtF);

			//���[�gG�����l����
			double min_G = ev_start2k_k2end(rt_G, 0, r_gene[gg][i]) + gr_k;
			int min_rt_G = rt_G[0];//�f�|�̊ԁi�f�|�ƍŏ��̓s�s�̊ԁj�ɑ}��
			bool insert_gr_startG = true;
			each_route_min(rt_G, i, gr_k, countG, min_G, min_rt_G, insert_gr_startG, insert_end_rtG);*/


			//A�EB�EC�̕]������ׂ�
			ev_sort[0] = min_A;
			ev_sort[1] = min_B;
			ev_sort[2] = min_C;
			/*ev_sort[3] = min_D;
			ev_sort[4] = min_E;
			ev_sort[5] = min_F;
			ev_sort[6] = min_G;*/

			evABC_sort(ev_sort);


			if (min_A == ev_sort[0]) {

				insert_position(rt_A,countA,insert_gr_startA, insert_end_rtA, r_gene[gg][i], min_rt_A);

			}
			else if (min_B == ev_sort[0]) {
				insert_position(rt_B, countB, insert_gr_startB, insert_end_rtB, r_gene[gg][i], min_rt_B);
				
			}
			else if(min_C == ev_sort[0]) {
				insert_position(rt_C, countC, insert_gr_startC, insert_end_rtC, r_gene[gg][i], min_rt_C);
				
			}
			/*else if (min_D == ev_sort[0]) {
				insert_position(rt_D, countD, insert_gr_startD, insert_end_rtD, r_gene[gg][i], min_rt_D);

			}
			else if (min_E == ev_sort[0]) {
				insert_position(rt_E, countE, insert_gr_startE, insert_end_rtE, r_gene[gg][i], min_rt_E);

			}
			else if (min_F == ev_sort[0]) {
				insert_position(rt_F, countF, insert_gr_startF, insert_end_rtF, r_gene[gg][i], min_rt_F);

			}
			else if (min_G == ev_sort[0]) {
				insert_position(rt_G, countG, insert_gr_startG, insert_end_rtG, r_gene[gg][i], min_rt_G);

			}*/

			insert_end_rtA = false;
			insert_end_rtB = false;
			insert_end_rtC = false;
			insert_end_rtD = false;
			insert_end_rtE = false;
			insert_end_rtF = false;
			insert_end_rtG = false;
			insert_gr_startA = false;
			insert_gr_startB = false;
			insert_gr_startC = false;
			/*insert_gr_startD = false;
			insert_gr_startE = false;
			insert_gr_startF = false;
			insert_gr_startG = false;*/
		}

		

		if (countA > 1) {
			//cout << "A�ˑR�ψ�" << endl;
			mutation(rt_A, countA);//�ˑR�ψ�
		}
		
		if (countB > 1){
			//cout << "B�ˑR�ψ�" << endl;
			mutation(rt_B, countB);//�ˑR�ψ�
		}
		
		if (countC > 1) {
			//cout << "C�ˑR�ψ�" << endl;
			mutation(rt_C, countC);//�ˑR�ψ�
		}

		/*if (countD > 1) {
			//cout << "C�ˑR�ψ�" << endl;
			mutation(rt_D, countD);//�ˑR�ψ�
		}
		if (countE > 1) {
			//cout << "C�ˑR�ψ�" << endl;
			mutation(rt_E, countE);//�ˑR�ψ�
		}
		if (countF > 1) {
			//cout << "C�ˑR�ψ�" << endl;
			mutation(rt_F, countF);//�ˑR�ψ�
		}
		if (countG > 1) {
			//cout << "C�ˑR�ψ�" << endl;
			mutation(rt_G, countG);//�ˑR�ψ�
		}*/
		
		//2�̌o�H���Ƃ�2opt�ߖT
		if (countA > 1 && countB > 1 && countC > 1) {
			two_route_search(rt_A, rt_B, rt_C, rt_D, rt_E, rt_F, rt_G);
		}

		cout << "9�E�E�E�E�E�E�E�E�E" << endl;

		//���̃R�X�g
		
		dist_A = dist_ABC(rt_A, countA);
		dist_B = dist_ABC(rt_B, countB);
		dist_C = dist_ABC(rt_C, countC);
		/*dist_D = dist_ABC(rt_D, countD);
		dist_E = dist_ABC(rt_E, countE);
		dist_F = dist_ABC(rt_F, countF);
		dist_G = dist_ABC(rt_G, countG);*/

		if (A == 3) {
			now_cost = cost3(dist_A, dist_B, dist_C);
		}
		else if (A == 5) {
			now_cost = cost5(dist_A, dist_B, dist_C,dist_D,dist_E);
		}
		else {
			now_cost = cost7(dist_A, dist_B, dist_C, dist_D, dist_E,dist_F,dist_G);
		}

		if (now_cost < min_cost) {
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
			for (int p = 0; p < countF + 1; p++) {
				best_F[p] = rt_F[p];
			}
			for (int p = 0; p < countG + 1; p++) {
				best_G[p] = rt_G[p];
			}*/

			best_countA = countA;
			best_countB = countB;
			best_countC = countC;
			/*best_countD = countD;
			best_countE = countE;
			best_countF = countF;
			best_countG = countG;*/

		}
		cout << gg << endl;
		gg++;


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

	//�f�|�̕\��
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

	/* �r�b�g�}�b�v������̕`�� */
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

double ev_start2k_k2end(int rt_ABC[], int start_or_end,int rt) {
	double gr_ABC, ABC_k;

	gr_ABC = sqrt(pow((gr_x - pos[rt_ABC[start_or_end]][0]), 2) + pow((gr_y - pos[rt_ABC[start_or_end]][1]), 2));

	ABC_k = sqrt(pow((pos[rt_ABC[start_or_end]][0] - pos[rt][0]), 2) + pow((pos[rt_ABC[start_or_end]][1] - pos[rt][1]), 2));

	return ABC_k - gr_ABC;

}

void render_string(float x, float y, const char* str, double cost) { //�����\�L�p
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

double cost3(double dist_A,double dist_B,double dist_C) {
	double sum = 0, ave = 0;

	sum = dist_A + dist_B + dist_C;
	ave = sum / A;

	return sum + aa*(fabs(ave - dist_A) + fabs(ave - dist_B) + fabs(ave - dist_C))/A;

}

double cost5(double dist_A, double dist_B, double dist_C, double dist_D, double dist_E) {
	double sum = 0, ave = 0;

	sum = dist_A + dist_B + dist_C + dist_D + dist_E;
	ave = sum / A;

	return sum + aa * (fabs(ave - dist_A) + fabs(ave - dist_B) + fabs(ave - dist_C)+ fabs(ave - dist_D)+ fabs(ave - dist_E)) / A;

}

double cost7(double dist_A, double dist_B, double dist_C, double dist_D, double dist_E, double dist_F, double dist_G) {
	double sum = 0, ave = 0,dif=0;

	sum = dist_A + dist_B + dist_C + dist_D + dist_E + dist_F + dist_G;
	ave = sum / A;
	dif = fabs(ave - dist_A) + fabs(ave - dist_B) + fabs(ave - dist_C) + fabs(ave - dist_D) + fabs(ave - dist_E) + fabs(ave - dist_F) + fabs(ave - dist_G);
	
	return sum + aa * dif / A;

}

double dist_ABC(int route[], int count) {
	double distance = 0;

	distance += sqrt(pow((gr_x - pos[route[0]][0]), 2) + pow((gr_y - pos[route[0]][1]), 2));
	for (int i = 0; i < count; i++) {
		distance += sqrt(pow((pos[route[i]][0] - pos[route[i + 1]][0]), 2) + pow((pos[route[i]][1] - pos[route[i + 1]][1]), 2));
	}
	distance += sqrt(pow((gr_x - pos[route[count]][0]), 2) + pow((gr_y - pos[route[count]][1]), 2));

	return distance;
}

//2�̃��[�g�ǂ����ł̋���
double dist_two_route(int route[], int count) {
	double distance = 0;	
	for (int i = 0; i < count; i++) {
		distance += sqrt(pow((pos[route[i]][0] - pos[route[i + 1]][0]), 2) + pow((pos[route[i]][1] - pos[route[i + 1]][1]), 2));
	}
	
	return distance;
}

void mutation(int *route,int count) {
	double dist_ac, dist_cb,dist_ab,now_dist;// , dist_pn, dist_ab, dist_pc, dist_cn, now_dist;
	double min_mutation = 100;//�ω��ʂ��ŏ��Ȃ��̂�ێ�
	double now_mutation = 0;//���̕ω���
	int insert_city = 0;//�}������s�s
	int pre_city = 0;//���̓s�s�̎��ɑ}��
	bool is_improved = true;//���P���ꂽ��
	bool is_improved_2opt = true;
	double min_dist = dist_ABC(route, count);//�Ǐ��T������O�̋���

	//1.5-opt�ߖT
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
			for (int c = k + 2; c < count+1; c++) {
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

			for (int s = insert_city; s > pre_city+1; s--) {
				swap(route[s], route[s - 1]);//�O�̔z��ƌ���
			}
		}
		else {
			for (int s = insert_city; s < pre_city; s++) {
				swap(route[s], route[s + 1]);//���̔z��ƌ���
			}
		}

		now_dist = dist_ABC(route, count);

		if (now_dist < min_dist) {
			min_dist = now_dist;
			//cout << "1.5-opt�ߖT���s" << endl;
			is_improved = true;
		}
		else {
			//���Ƃɖ߂�
			if (insert_city > pre_city) {
				for (int s = pre_city + 1; s < insert_city; s++) {
					swap(route[s], route[s + 1]);//�O�̔z��ƌ���
				}
			}
			else {
				for (int s = insert_city; s < pre_city; s++) {
					swap(route[s], route[s + 1]);//�O�̔z��ƌ���
				}
			}
			is_improved = false;
		}
	}
	
	//2-opt�ߖT
	int now_temp[N] = {0},min_x=0,min_c=0;
	//�l���R�s�[
	for (int x = 0; x < count + 1; x++) {
		now_temp[x] = route[x];
	}
	is_improved = true;
	min_dist = dist_ABC(route, count);
	double pre_dist = dist_ABC(route, count);//�Ǐ��T���@�i2-opt�ߖT�j�̑O�̋���
	
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
				//���ɖ߂�
				for (int x = 0; x < count + 1; x++) {
					now_temp[x] = route[x];
				}
			}
		}
		
		if (min_dist < pre_dist) {
			pre_dist = min_dist;
			//�o�H���X�V
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
	
	

}



void swap(int &x, int &y) {
	int temp;    // �l���ꎞ�ۑ�����ϐ�

	temp = x;
	x = y;
	y = temp;
}

int GetRandom(int min, int max) {

	return min + (int)(rand() * (max - min + 1.0) / (1.0 + RAND_MAX));

}

double insert_but_gr(int route[],int insert_city,int j) {
	double X_Y=0, X_k=0, Y_k=0;

	X_Y = sqrt(pow((pos[route[j]][0] - pos[route[j + 1]][0]), 2) + pow((pos[route[j]][1] - pos[route[j + 1]][1]), 2));
	X_k = sqrt(pow((pos[route[j]][0] - pos[r_gene[gg][insert_city]][0]), 2) + pow((pos[route[j]][1] - pos[r_gene[gg][insert_city]][1]), 2));
	Y_k = sqrt(pow((pos[route[j + 1]][0] - pos[r_gene[gg][insert_city]][0]), 2) + pow((pos[route[j + 1]][1] - pos[r_gene[gg][insert_city]][1]), 2));

	return X_k + Y_k - X_Y;
}

void each_route_min(int route[],int insert_city,double gr_k,int count,double& min, int &min_rt,bool &insert_gr_start,bool &insert_end_rt) {

	double ev_XY, ev_end;
	

	if (count != 0) {
		for (int j = 0; j < count; j++) {
			//���[�g���ő}���i�f�|�ȊO�j
			ev_XY = insert_but_gr(route, insert_city, j);

			if (min > ev_XY) {
				min = ev_XY;
				min_rt = j;//���̓s�s�̎��ɂ���}��
				insert_gr_start = false;

			}
		}

		//�f�|�ɋA���Ă���O�ɑ}��
		ev_end = ev_start2k_k2end(route, count, r_gene[gg][insert_city]) + gr_k;

		if (min > ev_end) {
			min = ev_end;

			insert_end_rt = true;
			insert_gr_start = false;

		}

	}


}

void insert_position(int route[],int &count,bool insert_gr_start,bool insert_end_rt,int insert_city,int min_rt) {
	//�f�|�ƍŏ��̓s�s�̊Ԃɑ}��
	if (insert_gr_start) {
		count++;
		routeInsert(route, count, insert_city, 0);

	}
	//�Ō�ƃf�|�̊Ԃɑ}��
	else if (insert_end_rt) {
		count++;
		route[count] = insert_city;

	}
	else {
		count++;
		routeInsert(route, count, insert_city, min_rt + 1);

	}
}

//2-opt���s��2�̃��[�g�̑g����
void two_route_search(int* rt_A, int* rt_B, int* rt_C, int* rt_D, int* rt_E, int* rt_F, int* rt_G ) {
	cout << "aaa�E�E�E�E�E�E�E�E�E" << endl;
	if (num_car == 3) {
		//att48[3]�̂Ƃ��Ɏ��s
		cout << "A�E�E�E�E�E�E�E�E�E" << endl;
		search_in_each_route(rt_A, rt_B, countA, countB);
		cout << "B�E�E�E�E�E�E�E�E�E" << endl;
		search_in_each_route(rt_B, rt_C, countB, countC);
		cout << "C�E�E�E�E�E�E�E�E�E" << endl;
		search_in_each_route(rt_C, rt_A, countC, countA);
		
	}
	else if (num_car == 5) {
		//pcb442[5]�̂Ƃ��Ɏ��s
		search_in_each_route(rt_A, rt_B, countA, countB);
		search_in_each_route(rt_B, rt_C, countB, countC);
		search_in_each_route(rt_C, rt_D, countC, countD);
		search_in_each_route(rt_D, rt_E, countD, countE);
		search_in_each_route(rt_E, rt_A, countE, countA);
	}
	else {
		//pr2392[7]�̂Ƃ��Ɏ��s
		search_in_each_route(rt_A, rt_B, countA, countB);
		search_in_each_route(rt_B, rt_C, countB, countC);
		search_in_each_route(rt_C, rt_D, countC, countD);
		search_in_each_route(rt_D, rt_E, countD, countE);
		search_in_each_route(rt_E, rt_F, countE, countF);
		search_in_each_route(rt_F, rt_G, countF, countG);
		search_in_each_route(rt_G, rt_A, countG, countA);
	}
	
}

void search_in_each_route(int* rt_1,int* rt_2,int count_1,int count_2) {
	int split_1 = count_1+1;		//���[�g�𕪂���Ƃ��p
	int split_2 = split_1+count_2+1;	//���[�g�𕪂���Ƃ��p
	int num = count_1+count_2;
	int newtemp[N] = { 0 };
	int start_2opt_position = 0;	//2opt������͈͂̃X�^�[�g
	int end_2opt_position = 0;		//2opt������͈͂̃G���h

	double min_dist = 0,now_dist=0;
	cout << "1�E�E�E�E�E�E�E�E�E" << endl;
	//3�̃��[�g����̃��[�g��
	for (int i = 0; i < split_1; i++) {
		newtemp[i] = rt_1[i];
	}
	for (int i = split_1; i < split_2; i++) {
		newtemp[i] = rt_2[i - split_1];
	}
	cout << "2�E�E�E�E�E�E�E�E�E" << endl;
	//���݂̋����i2�̃��[�g�̋����j���ŒZ�Ƃ���
	min_dist = dist_ABC(newtemp, num);
	cout << "3�E�E�E�E�E�E�E�E�E" << endl;
	
	//2opt�ߖT�̎��s
	for (int s = count_1 - 1; s <= count_1; s++) {
		for (int e = count_1 + 2; e > count_1; e--) {
			cout << "4�E�E�E�E�E�E�E�E�E" << endl;

			for (int i = s, j = e; i < count_1 + 1; i++, j--) {
				swap(newtemp[i], newtemp[j]);
			}
			
			//���s��̋��������߂�
			now_dist = dist_ABC(newtemp, num);

			//�ŒZ�ɂł��邩
			if (now_dist < min_dist) {
				//�ŒZ�ł���Ƃ���2opt������͈͂̂͂��߂Ƃ�����ێ�����
				start_2opt_position = s;
				end_2opt_position = e;
				cout << "���s�E�E�E�E�E�E�E�E�E" << endl;
			}

			//�t�������ɖ߂�
			for (int i = s, j = e; i < count_1 + 1; i++, j--) {
				swap(newtemp[i], newtemp[j]);
			}
			now_dist = 0;
		}

	}
	cout << "5�E�E�E�E�E�E�E�E�E" << endl;

	//���̃��[�g���ŒZ�̃��[�g�ɕύX
	for (int i = start_2opt_position, j = end_2opt_position; i < count_1 + 1; i++, j--) {
		swap(newtemp[i], newtemp[j]);
	}
	cout << "6�E�E�E�E�E�E�E�E�E" << endl;

	//�܂Ƃ߂����[�g��2�̃��[�g�ɖ߂��A�ŒZ�ɍX�V
	for (int i = 0; i < split_1; i++) {
		rt_1[i] = newtemp[i];
	}
	cout << "7�E�E�E�E�E�E�E�E�E" << endl;

	for (int j = split_1; j < split_2; j++) {
		rt_2[j - split_1] = newtemp[j];
	}
	cout << "8�E�E�E�E�E�E�E�E�E" << endl;

}