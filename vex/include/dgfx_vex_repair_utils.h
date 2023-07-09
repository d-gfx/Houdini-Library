/**
 * Copyright (c) 2023 d-gfx
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */
#ifndef dgfx_vex_repair_utils_h
#define dgfx_vex_repair_utils_h

#include "dgfx_vex_utils.h"

/**
 * Deprecated : Scheduled for deletion because it can be achieved by combining existing nodes
 * reference : https://dlanggraphics.blogspot.com/2023/07/HoudiniPolyDoctor001_01927828638.html
 * 
 * if Vertex is degeneracy, remove it
 * how to use : in Vertex Wrangle
 *  #include "dgfx_vex_utils.h"
 *  int geo = 0, is_dbg = 0;
 *  dgfx_Remove_Degeneracy_Vertex(geo, @vtxnum, @P, is_dbg);
 */
struct dgfx_Remove_Degeneracy_Vertex
{
    int is_dbg = 0;
    float eps = 0.0001;

    void remove(const int geo, vtxnum; const vector P)
    {
        int prim = vertexprim(geo, vtxnum);
        int nvtx = primvertexcount(geo, prim);
        // linear vertex number -> prim vertex index
        int self_ivtx = vertexprimindex(geo, vtxnum);

        int vtx_next = (self_ivtx + 1) % nvtx;
        int vtx_prev = (self_ivtx - 1) % nvtx;
        // convert linear vertex number
        vtx_next = vertexindex(geo, prim, vtx_next);
        vtx_prev = vertexindex(geo, prim, vtx_prev);

        vector vtx_next_P = vertex(geo, "P", vtx_next);
        vector vtx_prev_P = vertex(geo, "P", vtx_prev);
        if (is_dbg)
        {
            setvertexattrib(geo, "vtx_num", prim, self_ivtx, vtxnum);
            setvertexattrib(geo, "vtx_next", prim, self_ivtx, vtx_next);
            setvertexattrib(geo, "vtx_prev", prim, self_ivtx, vtx_prev);
        }

        if ((ptlined(vtx_prev_P, P, vtx_next_P) < eps)
        || (ptlined(vtx_next_P, P, vtx_prev_P) < eps))
        {
            if (is_dbg)
            {
                int pt = vertexpoint(geo, vtxnum);
                printf("remove vtx[%d] : pt[%d], prim[%d]\n", vtxnum, pt, prim);
            }
            removevertex(geo, vtxnum);
        }
    }
};

/**
 * Repair "crossed boundary (unshared) edges"
 *   that occur when using Labs Lot Subdivision, etc.
 * Define as a member function of a structure since it consists of 3 steps.
 *   1st : Primitive Wrangle, 2nd : Point Wrangle, 3rd : Primitive Wrangle
 * Temporary attribute names should be prefixed with "rcbue".
 */
struct dgfx_Repair_Crossed_Boundary_Unshared_Edges
{
    int is_dbg = 0;
    float eps = 0.001;
    // 1st step
    void Step_1st_Primitive_Wrangle(const int geo, primnum)
    {
        // Prepare an empty array of attributes to store.
        int append_pts[];
        setprimattrib(geo, "rcbue__append_pts", primnum, append_pts);
        int pts[] = primpoints(geo, primnum);
        if (is_dbg != 0)
        {
            setprimattrib(geo, "rcbue__orig_prim_pts", primnum, pts);
        }
    }
    // 2nd step
    void Step_2nd_Point_Wrangle(const int geo, ptnum; const vector P)
    {
        int pt_prims[] = pointprims(geo, ptnum);
        string prim_grp = "*";
        // Exclude primitives to which the Point already belongs from the search.
        foreach (int pt_prim; pt_prims)
        {
            prim_grp = sprintf("%s ^%s", prim_grp, pt_prim);
        }
        // Search for the nearest primitive to which the Point does not belong.
        int prim = -1; vector uvw;
        float dist = xyzdist(geo, prim_grp, P, prim, uvw);
        vector hit_P = primuv(geo, "P", prim, uvw);
        // Determines if the specified position P is on the edge of the primitive found by xyzdist.
        if (!dgfx_Pos_Is_On_Primitive_Edges(geo, prim, hit_P, eps))
        {
            prim = -1;
        }
        else
        {
            setprimattrib(geo, "rcbue__append_pts", prim, ptnum, "append");
        }
        if (is_dbg != 0)
        {
            setprimattrib(geo, "rcbue__dist", prim, dist, "set");
            setpointattrib(geo, "rcbue__find_prim", ptnum, prim);
            setpointattrib(geo, "rcbue__prim_grp", ptnum, prim_grp);
        }
    }
    // 3rd step
    void Step_3rd_Primitive_Wrangle(const int geo, primnum)
    {
        int pts[] = primpoints(geo, primnum);
        int append_pts[] = prim(geo, "rcbue__append_pts", primnum);
        if (len(append_pts) == 0) { return; }
        // Check which edge the Point to be added is on,
        //   and add it to the position of the array
        foreach (int append_pt; append_pts)
        {
            int npts = len(pts);
            vector append_P = point(geo, "P", append_pt);
            for (int i=0; i<npts; ++i)
            {
                vector pt1_P = point(geo, "P", pts[i]);
                vector pt2_P = point(geo, "P", pts[(i+1)%npts]);
                float dist = ptlined(pt1_P, pt2_P, append_P);
                if (dist < eps)
                {
                    insert(pts, i+1, append_pt);
                    break;
                }
            }
        }
        int add_prim = addprim(geo, "poly", pts);
        removeprim(geo, primnum, 0);
        if (is_dbg != 0)
        {
            setprimattrib(geo, "rcbue__new_prim_pts", add_prim, pts);
            setprimattrib(geo, "rcbue__orig_prim", add_prim, primnum);
        }
    }
};

#endif // dgfx_vex_repair_utils_h
