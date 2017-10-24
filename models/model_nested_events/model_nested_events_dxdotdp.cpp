
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_nested_events_dwdp.h"
#include "model_nested_events_w.h"

int dxdotdp_model_nested_events(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
int ip;
int ix;
memset(tdata->dxdotdp,0,sizeof(realtype)*1*udata->nplist);
status = dwdp_model_nested_events(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 3: {
  tdata->dxdotdp[0 + ip*model->nx] = tdata->h[1]*x_tmp[0];

  } break;

  case 4: {
  tdata->dxdotdp[0 + ip*model->nx] = -x_tmp[0];

  } break;

}
}
for(ip = 0; ip<udata->nplist; ip++) {
   for(ix = 0; ix<model->nx; ix++) {
       if(amiIsNaN(tdata->dxdotdp[ix+ip*model->nx])) {
           tdata->dxdotdp[ix+ip*model->nx] = 0;
           if(!tdata->nan_dxdotdp) {
               warnMsgIdAndTxt("AMICI:mex:fdxdotdp:NaN","AMICI replaced a NaN value in dxdotdp and replaced it by 0.0. This will not be reported again for this simulation run.");
               tdata->nan_dxdotdp = TRUE;
           }
       }
       if(amiIsInf(tdata->dxdotdp[ix+ip*model->nx])) {
           warnMsgIdAndTxt("AMICI:mex:fdxdotdp:Inf","AMICI encountered an Inf value in dxdotdp, aborting.");
           return(-1);
       }
   }
}
return(status);

}

