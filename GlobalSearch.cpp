/**
 * @file GlobalSearch.cpp
 * @author Wang Xi 王熙 (xiwang_chn@foxmail.com)
 * @brief It provides an implementation of the meta contact search （MCS) algorithm in article
 * ”Wang, X., Wu, W., Zhu, H., Zhang, H., Lin, J. S., & Bobet, A. (2022). A global direct search method for high-fidelity contact detection between arbitrarily shaped three-dimensional convex polyhedral blocks. Computers and Geotechnics, 150, 104891."
 * @version 0.1
 * @date 2023-10-25
 *
 * @copyright Copyright (c) 2023, see the LICENSE, README.md
 */

#include "GlobalSearch.h"

int DetectPolyCtBlkGlb(
    Package *pkg, Env *env, uint nB1, uint nB2, Contacts *contacts)
{
  // Get data
  // Get polyhedral index data
  Polyhedron *po1 = &pkg->polyhedra[nB1],
             *po2 = &pkg->polyhedra[nB2];
  // Get vertex array
  const float3h *vertexArr1 = GetBlockNodeArrConst(pkg, po1),
                *vertexArr2 = GetBlockNodeArrConst(pkg, po2);
  // Get map from vertex to edges
  const uchar *v2eMapArr1 = GetV2EMapArrConst(pkg, po1);
  const uchar *v2eMapArr2 = GetV2EMapArrConst(pkg, po2);
  // Get Edge array, it includes node indexes
  const Edge *edgeArr1 = GetEdgeArrConst(pkg, po1),
             *edgeArr2 = GetEdgeArrConst(pkg, po2);
  // Get Face array, it includes node indexes
  const Face *faceArr1 = GetFaceFeatureArrConst(pkg, po1),
             *faceArr2 = GetFaceFeatureArrConst(pkg, po2);

  // valid: strictly valid; loose, allow for some tolerance, see the article Wang et al.(2022)
  vecPolyCt contacts_valid, contacts_loose;
  contacts_valid.reserve(MAX(po1->nNodeLen, po2->nNodeLen));

  // loop to compute VF
  if (!DetVFGlb(nB1, nB2, po1->nNodeLen, po2->nFaceLen,
                vertexArr1, edgeArr1, v2eMapArr1,
                vertexArr2, faceArr2, pkg->faceNodeList,
                &pkg->boxs[nB2], contacts_valid, contacts_loose))

    return 0;
  // // loop to compute FV
  if (!DetVFGlb(nB2, nB1, po2->nNodeLen, po1->nFaceLen,
                vertexArr2, edgeArr2, v2eMapArr2,
                vertexArr1, faceArr1, pkg->faceNodeList,
                &pkg->boxs[nB1], contacts_valid, contacts_loose))
    return 0;

  // loop to compute CEE
  if (!DetCEEGlb(nB1, nB2, po1->nEdgeLen, po2->nEdgeLen,
                 vertexArr1, edgeArr1, faceArr1,
                 vertexArr2, edgeArr2, faceArr2,
                 contacts_valid, contacts_loose))
    return 0;

  // select contacts, transfer contacts and compute relative velocity
  return FilterContacts(pkg, env, contacts_valid, contacts_loose, contacts);
}

bool DetVFGlb(
    uint nB1, uint nB2, uint nVLen1, uint nFLen2,
    const float3h *vertexArr1, const Edge *edgeArr1, const uchar *v2eMapArr1,
    const float3h *vertexArr2, const Face *faceArr2, const uchar *faceIdArr,
    const BlkBox *box2, vecPolyCt &contacts_valid, vecPolyCt &contacts_loose)
{
  for (size_t i = 0; i < nVLen1; i++)
  {
    const float3h *v1 = &vertexArr1[i];
    const uchar *v2eMap1 = &v2eMapArr1[i * NVEFM];

    for (size_t j = 0; j < nFLen2; j++)
    {
      /// OVERLAP:
      const Face *f2 = &faceArr2[j];
      bool strictNoOverlap = true;    // to judge if no contact
      bool strictNoOverlap_Ct = true; // to select contact for direction
      uchar k;
      for (k = 1; k <= v2eMap1[0]; k++)
      {
        const Edge *E = &edgeArr1[v2eMap1[k]];
        const float3h &v12 = i == E->eNode1 ? E->info.e12V : -1 * E->info.e12V;
        // if do not satisfy the overlapping condition, skip this VF pair
        float_ dOverlap = Dot(v12, f2->normal);
        if (dOverlap < -dOverlapAng)
          break;
        if (dOverlap < 0)
          strictNoOverlap = false;
        if (dOverlap <= -SOVERA)
          strictNoOverlap_Ct = false;
      }
      // break early, fail the no-overlap check， skip this VF pair
      if (k <= v2eMap1[0])
        continue;

      /// DISTANCE:
      // no overlap, compute the distance, if >tol, return directly
      float_ dn = Dot(f2->normal, *v1 - f2->centroid);
      // to guarantee that no ghost contacts can exist when relaxing the overlap and penetration
      if (strictNoOverlap && dn >= 0)
        return false;

      /// POSITION:
      PolyhedraContact Ct;
      Ct.nB1 = nB1;
      Ct.nB2 = nB2;
      Ct.nCtType = 1;
      Ct.nP1 = i;
      Ct.nP2 = j;
      Ct.dn = dn;
      Ct.nDir = f2->normal;
      Ct.P1 = vertexArr1[i];
      Ct.P2 = Ct.P1 - Ct.dn * Ct.nDir;
      // Criterion of position
      Ct.bPos = NodeInConvexFace(*v1, vertexArr2,
                                 f2->length, &faceIdArr[f2->start]);

      /// SAVECT:
      if (strictNoOverlap_Ct)
        contacts_valid.push_back(Ct);
      else
        contacts_loose.push_back(Ct);
    }
  }
  return true;
}

// need to set a paralell
bool DetCEEGlb(
    uint nB1, uint nB2, uint nELen1, uint nELen2,
    const float3h *vertexArr1, const Edge *edgeArr1, const Face *faceArr1,
    const float3h *vertexArr2, const Edge *edgeArr2, const Face *faceArr2,
    vecPolyCt &contacts_valid, vecPolyCt &contacts_loose)
{
  for (size_t i = 0; i < nELen1; i++)
  {
    const float3h *e[4];
    const Edge *E1 = &edgeArr1[i];
    e[0] = &vertexArr1[E1->eNode1];
    e[1] = &vertexArr1[E1->eNode2];
    for (size_t j = 0; j < nELen2; j++)
    {
      const Edge *E2 = &edgeArr2[j];
      e[2] = &vertexArr2[E2->eNode1];
      e[3] = &vertexArr2[E2->eNode2];

      // four normals of the adjacent faces of two contacting edges
      const float3h *nf[4];
      nf[2] = &faceArr2[edgeArr2[j].nFace1].normal;
      nf[3] = &faceArr2[edgeArr2[j].nFace2].normal;
      // const float3h &v_1 = (*nf[0] + *nf[1]) / 2;
      const float3h &v_2 = (*nf[2] + *nf[3]) / 2;

      /// compute the right contact normal of CEE
      float3h n;
      Cross(E1->info.e12V, E2->info.e12V, n);
      float_ m_n = Modulus(n);
      if (m_n < PARAA) // if parallel, continue to skip numerical errors
        continue;

      /// NORMAL: OVERLAP: DIST:
      bool strictNoOverlap = true;
      bool strictNoOverlap_Ct = true; // to select contact for direction
      float_ dn;
      // Check normal
      n /= m_n;
      if (Dot(n, v_2) < 0)
        n = -1 * n;
      // check overlap
      float_ dOverlap10 = Dot(E1->info.I1, n), dOverlap20 = Dot(E1->info.I2, n),
             dOverlap31 = Dot(E2->info.I1, n), dOverlap41 = Dot(E2->info.I2, n);
      // do not satisfy no overlap
      if (dOverlap10 < -1 * dOverlapAng || dOverlap20 < -1 * dOverlapAng ||
          dOverlap31 > dOverlapAng || dOverlap41 > dOverlapAng)
        continue;
      // discern whether strict overlap or not
      strictNoOverlap = (dOverlap10 > 0 && dOverlap20 > 0 &&
                         dOverlap31 < 0 && dOverlap41 < 0);
      strictNoOverlap_Ct = (dOverlap10 > -SOVERA &&
                            dOverlap20 > -SOVERA &&
                            dOverlap31 < SOVERA &&
                            dOverlap41 < SOVERA);

      // dn<0 indicates penetration, dn>0 for open
      dn = Dot((*e[0] - *e[2]), n);

      // if strict no overlap, then the distance criterion should also be strict
      if (strictNoOverlap && dn > 0)
        return false;

      /// POSITION:
      float_ dU1 = 0, dU2 = 0;
      PolyhedraContact Ct;

      // get projection points on two contacting edges
      int nEE = EEPoints(e, dU1, dU2, Ct.P1, Ct.P2);
      Ct.bPos = nEE == 1;
      if (nEE < 0) // Consider as parallel
        continue;

      /// SAVECT:
      // store the contact index
      Ct.nB1 = nB1;
      Ct.nB2 = nB2;
      Ct.nCtType = 2;
      Ct.nP1 = i;
      Ct.nP2 = j;
      Ct.dn = dn;
      Ct.nDir = n;
      if (strictNoOverlap_Ct)
        contacts_valid.push_back(Ct);
      else
        contacts_loose.push_back(Ct);
    }
  }
  return true;
}

int FilterContacts(Package *pkg, Env *env,
                   vecPolyCt &contacts_valid, vecPolyCt &contacts_loose, Contacts *contacts)
{
  // select the contact with the max d
  // then select remaining contacts that shares the same contact normal
  int nMaxCt(-1);
  float_ dMaxD(-INFINITY);
  for (size_t i = 0; i < contacts_valid.size(); i++)
    if (contacts_valid[i].dn > dMaxD)
    {
      dMaxD = contacts_valid[i].dn;
      nMaxCt = i;
    }
  if (nMaxCt < 0) // detect some contact
  {
    // loose should also be absent
    assert(contacts_loose.size() == 0);
    return 0;
  }

  int nCt(0);
  uint nB1 = contacts_valid[nMaxCt].nB1;
  const float3h &n = contacts_valid[nMaxCt].nDir;

  for (size_t i = 0; i < contacts_valid.size(); i++)
  {
    PolyhedraContact &Ct = contacts_valid[i];
    if (!Ct.bPos) // Skip if the projection is not in the right place
      continue;
    if (Ct.dn < 0 && Ct.dn > -pkg->polyhedralSta.minInSphereD)
    {
      int sign = Ct.nB1 == nB1 ? 1 : -1;
      if (abs(Dot(Ct.nDir, sign * n) - 1) < dOverlapAng)
      {
        Ct.nMech = pkg->polyhedra[Ct.nB1].iMechFeature * pkg->nMechFeatures +
                   pkg->polyhedra[Ct.nB2].iMechFeature;
        // transfer contact status from the last step, compute relative velocity
        if (ContactTransfer(Ct, contacts) == 0) // no pre contact, use new
        {
          Ct.nCtSt = 2; // Initialize status as lock
          Ct.ds = 0;    // if first time, then ds=0
          Ct.sDir | 0;  // if first time, then sDir=0
        }
        const Polyhedron *po1 = &pkg->polyhedra[Ct.nB1];
        const Polyhedron *po2 = &pkg->polyhedra[Ct.nB2];
        ComputeRV(pkg, env, contacts, po1, po2, Ct);

#ifdef OPENMP
#pragma omp critical
#endif
        { // Save contacts, Atomic when parallel
          contacts->iPolyCt[Ct.nB1].push_back(contacts->vPolyCt.size());
          contacts->iPolyCt[Ct.nB2].push_back(contacts->vPolyCt.size());
          contacts->vPolyCt.push_back(Ct);
        }

        nCt++;
      }
    }
  }

  for (size_t i = 0; i < contacts_loose.size(); i++)
  {
    PolyhedraContact &Ct = contacts_loose[i];
    if (!Ct.bPos) // Skip if the projection is not in the right place
      continue;
    if (Ct.dn < 0 && Ct.dn > -pkg->polyhedralSta.minInSphereD)
    {
      int sign = Ct.nB1 == nB1 ? 1 : -1;
      if (abs(Dot(Ct.nDir, sign * n) - 1) < dOverlapAng)
      {
        Ct.nMech = pkg->polyhedra[Ct.nB1].iMechFeature * pkg->nMechFeatures +
                   pkg->polyhedra[Ct.nB2].iMechFeature;
        // transfer contact status from the last step, compute relative velocity
        if (ContactTransfer(Ct, contacts) == 0) // no pre contact, use new
        {
          Ct.nCtSt = 2; // Initialize status as lock
          Ct.ds = 0;    // if first time, then ds=0
          Ct.sDir | 0;  // if first time, then sDir=0
        }
        const Polyhedron *po1 = &pkg->polyhedra[Ct.nB1];
        const Polyhedron *po2 = &pkg->polyhedra[Ct.nB2];
        ComputeRV(pkg, env, contacts, po1, po2, Ct);

#ifdef OPENMP
#pragma omp critical
#endif
        { // Save contacts, Atomic when parallel
          contacts->iPolyCt[Ct.nB1].push_back(contacts->vPolyCt.size());
          contacts->iPolyCt[Ct.nB2].push_back(contacts->vPolyCt.size());
          contacts->vPolyCt.push_back(Ct);
        }

        nCt++;
      }
    }
  }

  return nCt; // return number of meta contacts between two convex polyhedral blocks
}
