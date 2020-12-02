#ifndef ROT_EULANG_H
#define ROT_EULANG_H

#include <stdbool.h>

/*! @file
* Euler angle ranges:                                                         
*                                                                             
* XYZ, ZYX: phi in [-pi, pi],     theta in [-pi/2, pi/2], psi in [-pi, pi]    
* XZY, YZX: phi in [-pi, pi],     theta in [-pi, pi],     psi in [-pi/2, pi/2]
* ZXY, YXZ: phi in [-pi/2, pi/2], theta in [-pi, pi],     psi in [-pi, pi]    
*                                                                             
* Euler angle sequence: XYZ (world). First rotation about X, second rotation  
* about Y, and the third rotation about Z axis of the world(i.e. fixed) frame.
* This is the same as the sequence used in Blender.                           
* In contrast, the XYZ sequence is understood in the Aerospace community as:  
* First rotation about Z-axis, second rotation about Y-axis, and the third    
* rotation about X-axis of the body frame.                                    
*                                                                             
* @see http://www.geometrictools.com/Documentation/EulerAngles.pdf            
*/

enum EulangSeq {ES_XYZ, ES_XZY, ES_YXZ, ES_YZX, ES_ZXY, ES_ZYX};

void rotmat_euler(const double* euler, enum EulangSeq seq, const bool world,
        double* rotmat);

void factor_rotmat(const double* rotmat, enum EulangSeq seq,
        const bool world, double* euler);


#endif //ROT_EULANG_H
