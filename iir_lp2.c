#define IMU_UPDATE_FREQ   500 
#define IIR_SHIFT         8

/**
 * Set ACC_WANTED_LPF1_CUTOFF_HZ to the wanted cut-off freq in Hz.
 * The highest cut-off freq that will have any affect is fs /(2*pi).
 * E.g. fs = 350 Hz -> highest cut-off = 350/(2*pi) = 55.7 Hz -> 55 Hz
 */
#define IMU_ACC_WANTED_LPF_CUTOFF_HZ  4

/**
 * Attenuation should be between 1 to 256.
 *
 * f0 = fs / 2*pi*attenuation ->
 * attenuation = fs / 2*pi*f0
 */
#define IMU_ACC_IIR_LPF_ATTENUATION (IMU_UPDATE_FREQ / (2 * 3.1415 * IMU_ACC_WANTED_LPF_CUTOFF_HZ))
#define IMU_ACC_IIR_LPF_ATT_FACTOR  (int)(((1<<IIR_SHIFT) / IMU_ACC_IIR_LPF_ATTENUATION) + 0.5)

typedef struct 
{
	s16 x;
	s16 y;
	s16 z;
} LPF;

typedef  struct
{
	s32 x;
	s32 y;
	s32 z;
} StoredFilterValues;

LPF accelLPF;
StoredFilterValues accelStoredFilterValues;

static u8 imuAccLpfAttFactor = IMU_ACC_IIR_LPF_ATT_FACTOR;

 /**
 * IIR filter the samples.
 */
int16_t iirLPFilterSingle(int32_t in, int32_t attenuation,  int32_t* filt)
{
  int32_t inScaled;
  int32_t filttmp = *filt;
  int16_t out;

  if (attenuation > (1<<IIR_SHIFT))
  {
    attenuation = (1<<IIR_SHIFT);
  }
  else if (attenuation < 1)
  {
    attenuation = 1;
  }

  // Shift to keep accuracy
  inScaled = in << IIR_SHIFT;
  // Calculate IIR filter
  filttmp = filttmp + (((inScaled-filttmp) >> IIR_SHIFT) * attenuation);
  // Scale and round
  out = (filttmp >> 8) + ((filttmp & (1 << (IIR_SHIFT - 1))) >> (IIR_SHIFT - 1));
  *filt = filttmp;

  return out;
}
                    
static void imuAccIIRLPFilter(xyz_s16_t* in, LPF* out, StoredFilterValues* storedValues, s32 attenuation)
{
  out->x = iirLPFilterSingle(in->x, attenuation, &storedValues->x);
  out->y = iirLPFilterSingle(in->y, attenuation, &storedValues->y);
  out->z = iirLPFilterSingle(in->z, attenuation, &storedValues->z);
}


imuAccIIRLPFilter(&, //int16  xyz
				  &accelLPF, //int16 xyz
				  &accelStoredFilterValues, //int32 storevalues xyz 
                 (int32_t)imuAccLpfAttFactor);// Attenuation  maybe error





#include <math.h>
/*
#define RATE_500_HZ 500
#define ATTITUDE_RATE RATE_500_HZ
#define ATTITUDE_RATE_LPF_CUTOFF_FREQ         30.0f
#define ATTITUDE_LPF_CUTOFF_FREQ              15.0f
*/
typedef struct {
  float a1;
  float a2;
  float b0;
  float b1;
  float b2;
  float delay_element_1;
  float delay_element_2;
} lpf2pData;


/**
 * 2-Pole low pass filter
 */
void lpf2pInit(lpf2pData* lpfData, float sample_freq, float cutoff_freq)
{
  if (lpfData == NULL || cutoff_freq <= 0.0f) {
    return;
  }

  lpf2pSetCutoffFreq(lpfData, sample_freq, cutoff_freq);
}

void lpf2pSetCutoffFreq(lpf2pData* lpfData, float sample_freq, float cutoff_freq)
{
  float fr = sample_freq/cutoff_freq;
  float ohm = tanf(M_PI_F/fr);
  float c = 1.0f+2.0f*cosf(M_PI_F/4.0f)*ohm+ohm*ohm;
  lpfData->b0 = ohm*ohm/c;
  lpfData->b1 = 2.0f*lpfData->b0;
  lpfData->b2 = lpfData->b0;
  lpfData->a1 = 2.0f*(ohm*ohm-1.0f)/c;
  lpfData->a2 = (1.0f-2.0f*cosf(M_PI_F/4.0f)*ohm+ohm*ohm)/c;
  lpfData->delay_element_1 = 0.0f;
  lpfData->delay_element_2 = 0.0f;
}

float lpf2pApply(lpf2pData* lpfData, float sample)//sample is acc
{
  float delay_element_0 = sample - lpfData->delay_element_1 * lpfData->a1 - lpfData->delay_element_2 * lpfData->a2;
  if (!isfinite(delay_element_0)) {
    // don't allow bad values to propigate via the filter
    delay_element_0 = sample;
  }

  float output = delay_element_0 * lpfData->b0 + lpfData->delay_element_1 * lpfData->b1 + lpfData->delay_element_2 * lpfData->b2;

  lpfData->delay_element_2 = lpfData->delay_element_1;
  lpfData->delay_element_1 = delay_element_0;
  return output;
}

float lpf2pReset(lpf2pData* lpfData, float sample)
{
  float dval = sample / (lpfData->b0 + lpfData->b1 + lpfData->b2);
  lpfData->delay_element_1 = dval;
  lpfData->delay_element_2 = dval;
  return lpf2pApply(lpfData, sample);
}
