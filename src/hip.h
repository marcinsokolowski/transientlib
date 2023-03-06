#define HIP_REC_LEN 58
#define HIPFIL_VJ 0
#define HIPFIL_VT 1
#define HIPFIL_BT 2
#define HIPFIL_VI 3
#define HIPFIL_BV 4
#define HIPFIL_I 5
 
struct hip_rec{
  double ra;
  double dec;
  long id;
  float vj,vt,vt_err,bt,bt_err,bv,bv_err,vi,vi_err;
  char src[2];
};

struct HIP{
  double ra,dec,mag;
  char src;
};
    

int SelectHipparcosStars(double alpha,double delta,
								 double equinox,double xbox,double ybox,
								 struct HIP **list,char* starcat,int filter);
