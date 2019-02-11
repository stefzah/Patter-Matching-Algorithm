#include < stdio.h > 
#include < stdlib.h >
  typedef struct {
    unsigned char r, g, b;
  }
pixel;
typedef struct {
  int x1, x2, y1, y2;
  double coef;
  int cul;
}
fereastra;
pixel ** CitireIMG(char * fisier) {
  int i, W, H, W2, H2, j;
  pixel ** P;
  unsigned char x;
  FILE * fin, * fout;
  fin = fopen(fisier, "rb");
  fseek(fin, 18, SEEK_SET);
  fread( & W, 4, 1, fin);
  fread( & H, 4, 1, fin);
  P = malloc(sizeof(int) * H);
  for (i = 0; i < H; i++)
    P[i] = malloc(sizeof(pixel) * W);
  fseek(fin, 54, SEEK_SET);
  for (i = H - 1; i >= 0; i--) {
    for (j = 0; j < W; j++) {
      fread( & x, 1, 1, fin);
      P[i][j].b = x;
      fread( & x, 1, 1, fin);
      P[i][j].g = x;
      fread( & x, 1, 1, fin);
      P[i][j].r = x;
    }
    j = W * 3;
    while (j % 4 != 0) {
      j++;
      fread( & x, 1, 1, fin);
    }
  }
  fclose(fin);
  return P;
}
double Calcul(int H, int W, int Hs, int Ws, fereastra f, pixel ** P, pixel ** Ps) {
  double Sp = 0, fp = 0, n = Ws * Hs, qS = 0, qf1 = 0, coeff = 0;
  int i, j;
  for (i = 0; i < Hs; i++)
    for (j = 0; j < Ws; j++)
      Sp += Ps[i][j].r;
  Sp /= Ws * Hs;
  for (i = f.x1; i <= f.x2; i++)
    for (j = f.y1; j <= f.y2; j++) {
      if (i >= 0 && j >= 0 && i < H && j < W) fp += P[i][j].r;
      else fp += 0;
    }
  fp /= Ws * Hs;
  for (i = 0; i < Hs; i++)
    for (j = 0; j < Ws; j++)
      qS += ((double) Ps[i][j].r - Sp) * ((double) Ps[i][j].r - Sp);
  qS /= (n + 1);
  qS = sqrt(qS);
  for (i = f.x1; i <= f.x2; i++)
    for (j = f.y1; j <= f.y2; j++) {
      if (i >= 0 && j >= 0 && i < H && j < W) qf1 += ((double) P[i][j].r - fp) * ((double) P[i][j].r - fp);
      else qf1 += fp * fp;
    }
  qf1 /= (n + 1);
  qf1 = sqrt(qf1);
  for (i = f.x1; i <= f.x2; i++)
    for (j = f.y1; j <= f.y2; j++) {
      if (i >= 0 && j >= 0 && i < H && j < W) coeff += 1 / (qf1 * qS) * ((double) P[i][j].r - fp) * ((double) Ps[i - f.x1][j - f.y1].r - Sp);
      else coeff += 1 / (qf1 * qS) * (-fp) * ((double) Ps[i - f.x1][j - f.y1].r - Sp);
    }
  coeff /= n;
  return coeff;
}
void Contur(char * fisier, fereastra f, pixel cul) {
  FILE * fin;
  int i, j, W, H;
  char c;
  pixel ** P, x;
  x.r = x.g = x.b = 0;
  fin = fopen(fisier, "r+b");
  fseek(fin, 18, SEEK_SET);
  fread( & W, 4, 1, fin);
  fread( & H, 4, 1, fin);
  P = CitireIMG(fisier);
  for (i = f.x1; i <= f.x2; i++) {
    j = f.y1;
    if (i >= 0 && j >= 0 && i < H && j < W) P[i][j] = cul;
    j = f.y2;
    if (i >= 0 && j >= 0 && i < H && j < W) P[i][j] = cul;
  }
  for (j = f.y1; j <= f.y2; j++) {
    i = f.x1;
    if (i >= 0 && j >= 0 && i < H && j < W) P[i][j] = cul;
    i = f.x2;
    if (i >= 0 && j >= 0 && i < H && j < W) P[i][j] = cul;
  }
  fseek(fin, 54, SEEK_SET);
  for (i = H - 1; i >= 0; i--) {
    for (j = 0; j < W; j++) {
      fwrite( & P[i][j].b, 1, 1, fin);
      fwrite( & P[i][j].g, 1, 1, fin);
      fwrite( & P[i][j].r, 1, 1, fin);
    }
    j = W * 3;
    while (j % 4 != 0) {
      j++;
      c = 0;
      fwrite( & c, 1, 1, fin);
    }
  }
  fclose(fin);
}
int Cifra(char * fisier) {
  if (strcmp(fisier, "cifra0.bmp") == 0) return 0;
  if (strcmp(fisier, "cifra1.bmp") == 0) return 1;
  if (strcmp(fisier, "cifra2.bmp") == 0) return 2;
  if (strcmp(fisier, "cifra3.bmp") == 0) return 3;
  if (strcmp(fisier, "cifra4.bmp") == 0) return 4;
  if (strcmp(fisier, "cifra5.bmp") == 0) return 5;
  if (strcmp(fisier, "cifra6.bmp") == 0) return 6;
  if (strcmp(fisier, "cifra7.bmp") == 0) return 7;
  if (strcmp(fisier, "cifra8.bmp") == 0) return 8;
  if (strcmp(fisier, "cifra9.bmp") == 0) return 9;
}
int * Coeficient(char * fisier, char * sablon, double prag) {
  int W, H, Ws, Hs, i, j, w, h, ii = 0, cif;
  fereastra f, * D;
  double calc;
  pixel ** P, ** Ps;
  unsigned char x;
  FILE * fin1, * fin2;
  fin1 = fopen(fisier, "rb");
  fin2 = fopen(sablon, "rb");
  cif = Cifra(sablon);
  fseek(fin1, 18, SEEK_SET);
  fread( & W, 4, 1, fin1);
  fread( & H, 4, 1, fin1);
  fseek(fin2, 18, SEEK_SET);
  fread( & Ws, 4, 1, fin2);
  fread( & Hs, 4, 1, fin2);
  P = CitireIMG(fisier);
  Ps = CitireIMG(sablon);
  D = malloc(sizeof(fereastra) * H * W);
  for (i = 0; i < H; i++)
    for (j = 0; j < W; j++) {
      f.cul = cif;
      w = Ws / 2;
      h = Hs / 2;
      if (Ws % 2 == 0) {
        f.y1 = j - w;
        f.y2 = j + w - 1;
      } else {
        f.y1 = j - w;
        f.y2 = j + w;
      }
      if (Hs % 2 == 0) {
        f.x1 = i - h;
        f.x2 = i + h - 1;
      } else {
        f.x1 = i - h;
        f.x2 = i + h;
      }
      calc = Calcul(H, W, Hs, Ws, f, P, Ps);
      if (prag <= calc) {
        ii++;
        D[ii] = f;
        D[ii].coef = calc;
      }
    }
  D = realloc(D, sizeof(fereastra) * (ii + 1));
  //printf("(%d)",ii+1);
  D[0].x1 = ii + 1;
  return D;
}
int cmp(const void * A,
  const void * B) {
  if (((fereastra * ) A) - > coef > ((fereastra * ) B) - > coef) return -1;
  else return 1;
}
void Sortare(fereastra * D) {
  qsort(D, D[0].x1, sizeof(fereastra), cmp);
}
double Arie(fereastra A) {
  return (double)(A.x2 - A.x1 + 1) * (A.y2 - A.y1 + 1);
}
double Intersectie(fereastra A, fereastra B) {
  double mn1, mn2, Ax1, Ax2, Bx1, Bx2, Ay1, By1, Ay2, By2;
  Ax1 = (double) A.x1;
  Ax2 = (double) A.x2;
  Ay1 = (double) A.y1;
  Ay2 = (double) A.y2;
  Bx1 = (double) B.x1;
  Bx2 = (double) B.x2;
  By1 = (double) B.y1;
  By2 = (double) B.y2;
  if (Ax2 < Bx2) mn1 = Ax2;
  else mn1 = Bx2;
  if (Ay2 < By2) mn2 = Ay2;
  else mn2 = By2;
  if (Ax1 <= Bx1 && Bx1 <= Ax2 && Ay1 <= By1 && By1 <= Ay2) return ((mn1 - Bx1 + 1) * (mn2 - By1 + 1));
  else if (Bx1 <= Ax1 && Ax1 <= Bx2 && By1 <= Ay1 && Ay1 <= By2) return ((mn1 - Ax1 + 1) * (mn2 - Ay1 + 1));
  else if (Bx1 <= Ax1 && Ax1 <= Bx2 && Ay1 <= By1 && By1 <= Ay2) return ((mn1 - Ax1 + 1) * (mn2 - By1 + 1));
  else if (Ax1 <= Bx1 && Bx1 <= Ax2 && By1 <= Ay1 && Ay1 <= By2) return ((mn1 - Bx1 + 1) * (mn2 - Ay1 + 1));
  else return 0;
}
void EliminareNonMaxime(fereastra * D) {
  int i, j, mx;
  for (i = 1; i < D[0].x1; i++) {
    if (D[i].coef != 0)
      for (j = i + 1; j < D[0].x1; j++) {
        if (D[j].coef != 0 && Intersectie(D[i], D[j]) / (Arie(D[i]) + Arie(D[j]) - Intersectie(D[i], D[j])) > 0.2) {
          D[j].coef = 0;
        }
      }
  }
  j = 1;
  mx = 0;
  for (i = 1; i < D[0].x1; i++, j++) {
    while (D[j].coef == 0) j++;
    if (j < D[0].x1) {
      D[i] = D[j];
      mx++;
    } else break;
  }
  D[0].x1 = mx + 1;
  D = realloc(D, sizeof(fereastra) * (mx + 1));
}
void AtribuireCoeficienti(char * fisier, fereastra ** D, int * ct, double coeff, int CT) {
  int i;
  char ** s;
  s = malloc(sizeof(int) * CT);
  for (i = 0; i < CT; i++) {
    s[i] = malloc(100);
    printf("Sablonul nr.%d: ", i);
    scanf("%s", s[i]);
    s[i] = realloc(s[i], strlen(s[i]));
  }
  for (i = 0; i < CT; i++) {
    D[i] = Coeficient(fisier, s[i], coeff);
    ( * ct) += D[i][0].x1 - 1;
  }
}
void AtribuireCulori(pixel * C) {
  C[0].r = 255;
  C[0].g = 0;
  C[0].b = 0;
  C[1].r = 255;
  C[1].g = 255;
  C[1].b = 0;
  C[2].r = 0;
  C[2].g = 255;
  C[2].b = 0;
  C[3].r = 0;
  C[3].g = 255;
  C[3].b = 255;
  C[4].r = 255;
  C[4].g = 0;
  C[4].b = 255;
  C[5].r = 0;
  C[5].g = 0;
  C[5].b = 255;
  C[6].r = 192;
  C[6].g = 192;
  C[6].b = 192;
  C[7].r = 255;
  C[7].g = 140;
  C[7].b = 0;
  C[8].r = 128;
  C[8].g = 0;
  C[8].b = 128;
  C[9].r = 128;
  C[9].g = 0;
  C[9].b = 0;
}
void PatternMatching(char * fisier) {
  fereastra ** D, * Dl;
  int i, j, * ct, CT = 0, CT1 = 0, NRS;
  pixel * C;
  ct = & CT;
  printf("Introduceti numarul de sabloane: ");
  C = malloc(sizeof(pixel) * 10);
  AtribuireCulori(C);
  scanf("%d", & NRS);
  D = malloc(NRS * sizeof(fereastra));
  AtribuireCoeficienti(fisier, D, ct, 0.5, NRS);
  CT++;
  Dl = malloc(sizeof(fereastra) * CT);
  CT = 1;
  for (i = 0; i < NRS; i++) {
    for (j = 1; j < D[i][0].x1; j++) {
      Dl[CT] = D[i][j];
      CT++;
    }
  }
  Dl[0].x1 = CT;
  Dl[0].coef = 1;
  Sortare(Dl);
  EliminareNonMaxime(Dl);
  CT = Dl[0].x1;
  for (i = 1; i < CT; i++) {
    Contur(fisier, Dl[i], C[Dl[i].cul]);
  }
  free(Dl);
  for (i = 0; i < NRS; i++)
    free(D[i]);
  free(D);
}
int main() {
  char * s;
  printf("Introduceti numele fisierului: ");
  s = malloc(100);
  scanf("%s", s);
  PatternMatching(s);
  free(s);
}
