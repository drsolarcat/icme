
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>
#include <vector>
#include <cstring>
#include "mex.h"

using namespace std;

tm* initTm() {
  time_t rawTime;
  ::time(&rawTime);
  return ::localtime(&rawTime);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  char *path;
  time_t beginUnixTime, endUnixTime, currentUnixTime;
  struct tm *currentTime;
  ifstream data;
  string line;
  istringstream linestr;
  vector<double> year, month, day, hour, minute, second, B, Bx, By, Bz,
    Vp, Vx, Vy, Vz, Pth, Np, Tp, Vth, beta;
  double tmpYear, tmpMonth, tmpDay, tmpHour, tmpMinute, tmpSecond,
    tmpB, tmpBx, tmpBy, tmpBz, tmpVp, tmpVx, tmpVy, tmpVz, tmpPth, tmpNp,
    tmpTp, tmpVth, tmpBeta;
  int n;

  path = mxArrayToString(prhs[0]);
  beginUnixTime = (mxGetScalar(prhs[1])-719529)*86400;
  endUnixTime = (mxGetScalar(prhs[2])-719529)*86400;
  currentTime = initTm();

  data.open(path, ifstream::in);
  if (data.is_open()) {
    while(data.good()) {
      getline(data, line);
      if (!line.empty() && line[0] != '#') {
        linestr.clear();
        linestr.str(line);
        linestr >> tmpYear >> tmpMonth >> tmpDay >> tmpHour >> tmpMinute >>
          tmpSecond >> tmpB >> tmpBx >> tmpBy >> tmpBz >> tmpVp >> tmpVx >>
          tmpVy >> tmpVz >> tmpPth >> tmpNp >> tmpTp >> tmpVth >> tmpBeta;
        currentTime->tm_year = tmpYear-1900;
        currentTime->tm_mon = tmpMonth-1;
        currentTime->tm_mday = tmpDay;
        currentTime->tm_hour = tmpHour;
        currentTime->tm_min = tmpMinute;
        currentTime->tm_sec = tmpSecond;
        currentUnixTime = mktime(currentTime);
        if (currentUnixTime > endUnixTime) break;
        if (currentUnixTime >= beginUnixTime &&
          currentUnixTime <= endUnixTime)
        {
          year.push_back(tmpYear);
          month.push_back(tmpMonth);
          day.push_back(tmpDay);
          hour.push_back(tmpHour);
          minute.push_back(tmpMinute);
          second.push_back(tmpSecond);
          B.push_back(tmpB);
          Bx.push_back(tmpBx);
          By.push_back(tmpBy);
          Bz.push_back(tmpBz);
          Vp.push_back(tmpVp);
          Vx.push_back(tmpVx);
          Vy.push_back(tmpVy);
          Vz.push_back(tmpVz);
          Pth.push_back(tmpPth);
          Np.push_back(tmpNp);
          Tp.push_back(tmpTp);
          Vth.push_back(tmpVth);
          beta.push_back(tmpBeta);
        }
      }
    }
    data.close();
  }

  n = year.size();

  double *pointers[] =
  {
    &year[0], &month[0], &day[0], &hour[0], &minute[0], &second[0],
    &B[0], &Bx[0], &By[0], &Bz[0], &Vp[0], &Vx[0], &Vy[0], &Vz[0],
    &Pth[0], &Np[0], &Tp[0], &Vth[0], &beta[0]
  };

  for (int i=0; i<=18; i++) {
    plhs[i] = mxCreateDoubleMatrix(1, n, mxREAL);
    memcpy(mxGetPr(plhs[i]), pointers[i], sizeof(double)*n);
  }
}

