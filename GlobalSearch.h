/**
 * @file GlobalSearch.cpp
 * @author Wang Xi 王熙 (xiwang_chn@foxmail.com)
 * @brief It provides an implementation of the algorithm in article:
 * ”Wang, X., Wu, W., Zhu, H., Zhang, H., Lin, J. S., & Bobet, A. (2022). A global direct search method for high-fidelity contact detection between arbitrarily shaped three-dimensional convex polyhedral blocks. Computers and Geotechnics, 150, 104891."
 * @version 0.1
 * @date 2023-10-25
 *
 * @copyright Copyright (c) 2023, see the LICENSE, README.md
 */

#ifndef _GLOBALSEARCH_H
#define _GLOBALSEARCH_H

// the invoke function of the meta contact search （MCS) algorithm
int DetectPolyCtBlkGlb(
    Package *pkg, Env *env, uint nB1, uint nB2, Contacts *contacts);

// detec vertex-face meta contacts
bool DetVFGlb(
    uint nB1, uint nB2, uint nVLen1, uint nFLen2,
    const float3h *vertexArr1, const Edge *edgeArr1, const uchar *v2eMapArr1,
    const float3h *vertexArr2, const Face *faceArr2, const uchar *faceIdArr,
    const BlkBox *box2, vecPolyCt &contacts_valid, vecPolyCt &contacts_loose);

// detect cross edge-edge meta contacts
bool DetCEEGlb(
    uint nB1, uint nB2, uint nELen1, uint nELen2,
    const float3h *vertexArr1, const Edge *edgeArr1, const Face *faceArr1,
    const float3h *vertexArr2, const Edge *edgeArr2, const Face *faceArr2,
    vecPolyCt &contacts_valid, vecPolyCt &contacts_loose);

// filter all meta contacts to get the actual contact plane
int FilterContacts(Package *pkg, Env *env,
                   vecPolyCt &contacts_valid, vecPolyCt &contacts_loose, Contacts *contacts);

#endif
