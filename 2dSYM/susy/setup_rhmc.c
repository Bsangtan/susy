// -----------------------------------------------------------------
// For now simply copy in output
// Wrap it into a routine that can be swapped out in the future
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Coefficients for (Mdag M)^(-1 / 4) and (Mdag M)^(1 / 8),
// with error < 2e-5 in range [1e-7, 1e3]
// Note the relative sign between fractional powers!
void setup_rhmc() {
  Norder = DEGREE;
  node0_printf("RHMC Norder %d\n", Norder);

  //awk '/res_MD/{print("amp4["$2"]=",$3";")}' < out_15.d
  ampdeg4 = 9.5797060554725838e-02;
  amp4[0] = 1.7701746700099842e-06;
  amp4[1] = 5.8705983656937455e-06;
  amp4[2] = 1.9961158693570120e-05;
  amp4[3] = 6.9125367600088173e-05;
  amp4[4] = 2.4032965323696816e-04;
  amp4[5] = 8.3620125835371663e-04;
  amp4[6] = 2.9099006745502945e-03;
  amp4[7] = 1.0126504714418652e-02;
  amp4[8] = 3.5241454044660878e-02;
  amp4[9] = 1.2266034741624667e-01;
  amp4[10] = 4.2721681852328125e-01;
  amp4[11] = 1.4932820692676758e+00;
  amp4[12] = 5.3188766358452595e+00;
  amp4[13] = 2.0944763089672641e+01;
  amp4[14] = 1.4525770103354523e+02;

  // awk '/pole_MD/{print("shift4["$2"]=",$3";")}' < out_15.d
  shift4[0] = 3.1085594175442315e-08;
  shift4[1] = 3.2994455960441383e-07;
  shift4[2] = 1.9424842756552213e-06;
  shift4[3] = 1.0453359626231250e-05;
  shift4[4] = 5.5337819905761986e-05;
  shift4[5] = 2.9204178440857227e-04;
  shift4[6] = 1.5403300046437174e-03;
  shift4[7] = 8.1233558140562465e-03;
  shift4[8] = 4.2840454273820550e-02;
  shift4[9] = 2.2594500626442715e-01;
  shift4[10] = 1.1921171782283737e+00;
  shift4[11] = 6.3026182343759860e+00;
  shift4[12] = 3.3683411978650057e+01;
  shift4[13] = 1.9083658214156412e+02;
  shift4[14] = 1.5386784635765257e+03;

  //awk '/res_GR/{print("amp8["$2"]=",$3";")}' < out_15.d
  ampdeg8 = 3.2148873149863206e+00;
  amp8[0] = -2.2977600408751347e-09;
  amp8[1] = -1.6898103706901084e-08;
  amp8[2] = -1.1099658368596436e-07;
  amp8[3] = -7.2162146587729939e-07;
  amp8[4] = -4.6841070484595924e-06;
  amp8[5] = -3.0396303865820389e-05;
  amp8[6] = -1.9723870959636086e-04;
  amp8[7] = -1.2798599250624023e-03;
  amp8[8] = -8.3051856063983548e-03;
  amp8[9] = -5.3904877281192094e-02;
  amp8[10] = -3.5026088217184553e-01;
  amp8[11] = -2.2893521967679966e+00;
  amp8[12] = -1.5436668340425719e+01;
  amp8[13] = -1.2297861076048798e+02;
  amp8[14] = -2.6252652966414048e+03;

  //awk '/pole_GR/{print("shift8["$2"]=",$3";")}' < out_15.d
  shift8[0] = 5.5367335615411457e-08;
  shift8[1] = 4.6910257304582898e-07;
  shift8[2] = 2.6768223190551614e-06;
  shift8[3] = 1.4319657256375662e-05;
  shift8[4] = 7.5694473187855338e-05;
  shift8[5] = 3.9922490005559548e-04;
  shift8[6] = 2.1046795395127538e-03;
  shift8[7] = 1.1094832053548640e-02;
  shift8[8] = 5.8486687698920667e-02;
  shift8[9] = 3.0834388405073770e-01;
  shift8[10] = 1.6264534005778293e+00;
  shift8[11] = 8.6030459456576764e+00;
  shift8[12] = 4.6179583183155444e+01;
  shift8[13] = 2.6854965277696181e+02;
  shift8[14] = 2.6004158696112045e+03;

  // Test by zeroing out all amp4 and amp8, optionally setting Norder to 1
//  int i;
//  Norder = 1;
//  for (i = 0; i < Norder; i++) {
//    amp4[i] = 0;
//    amp8[i] = 0;
//  }
//  node0_printf("TEST VERSION: internal Norder %d\n", Norder);
//  ampdeg4 = 0;
//  amp4[0] = 1;
//  shift4[0] = 1e-3;
//  ampdeg8 = 1;
//  amp8[0] = 0;
//  shift8[0] = 1;

#ifdef DEBUG_CHECK
  int i;
  node0_printf("RHMC ampdeg4 %e\n", ampdeg4);
  for (i = 0; i < Norder; i++)
    node0_printf("RHMC params %d amp4 %e shift4 %e\n", i, amp4[i], shift4[i]);

  node0_printf("RHMC ampdeg8 %e\n", ampdeg8);
  for (i = 0; i < Norder; i++)
    node0_printf("RHMC params %d amp8 %e shift8 %e\n", i, amp8[i], shift8[i]);
#endif
}
// -----------------------------------------------------------------
