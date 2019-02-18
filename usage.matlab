---------------------------------------------------------------------------

AUTHOR

Seungryong Kim

---------------------------------------------------------------------------

Version 2.0 (14 Apr. 2016)

---------------------------------------------------------------------------

CONTACT

web   : http://diml.yonsei.ac.kr/~srkim/DASC/
email : srkim89@yonsei.ac.kr

---------------------------------------------------------------------------

* The code is provided for academic use only. Use of the code in any commercial or industrial related activities is prohibited. 
* If you use this code in your research please give a reference to

@InProceedings{Kim2015,
author = {Seungryong Kim and Dongbo Min and Bumsub Ham and Seungchul Ryu and Minh N. Do and Kwanghoon Sohn},
title = {DASC: Dense Adaptive Self-Correlation Descriptor for Multi-modal and Multi-spectral Correspondence},
booktitle = {Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR), IEEE},
year = {2015}
}

---------------------------------------------------------------------------

* Usage
1) mex mexDASC.cpp
2) setup SIFTflow code [2]
3) start main.m

---------------------------------------------------------------------------

* Parameters
M_half: half size of large window M
N_half: half size of large window N 
epsil: epsilon for FastGuidedFilter [3] 
downSize: downsize factor s for FastGuidedFilter [3]
sigma_s: for recursive filter (RF) [4]
sigma_r: for recursive filter (RF) [4]
iter: for recursive filter (RF) [4]

---------------------------------------------------------------------------

* Input and Output
Input: input image 1                (e.g., 'img1.png')
       input image 2                (e.g., 'img2.png')
Output: warped image from image 2	(e.g., 'warp2.png')
        flow result             	(e.g., 'flow.png')

---------------------------------------------------------------------------

[1] S. Kim, D. Min, B. Ham, S. Ryu, M. N. Do., and K. Sohn, DASC: Dense 
Adaptive Self-Correlation Descriptor for Multi-modal and Multi-spectral 
Correspondence, In Proc. of CVPR, 2015.

[2] C. Liu, J. Yuen, and A. Torralba. Sift flow: Dense correspondence
across scenes and its applications. IEEE TPAMI, 33(5), pp. 815-830, 2011.

[3] K. He and J. Sun, Fast Guided Filter, arXiv, 2015.

[4] S. L. Eduardo and M. M. Oliveira, Domain transform for edge-aware image 
and video processing, ACM ToG, 30(4), 2011.

---------------------------------------------------------------------------