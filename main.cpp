#include <stdio.h>
#include <stdlib.h>
#include "glpk.h"
#include <vector>
#include <iostream>
#include <cstdio>
#include <iomanip>
#include <cassert>
#include <set>

using namespace std;
double eps = 1e-4;

glp_prob *mip;
int global_cnt_vert;
vector<int> *ref_left;
vector<int> *ref_right;

void dfs(int vert, set<int> &component, vector<vector<int>> &edj, vector<char> &visit) {
    if (visit[vert]) {
        return;
    }
    visit[vert] = 1;
    component.insert(vert);
    for (int next: edj[vert]) {
        dfs(next, component, edj, visit);
    }
}

void find_cut(glp_tree *tree, void *info) {
    if (glp_ios_reason(tree) == GLP_IROWGEN) {
        glp_prob *prob = glp_ios_get_prob(tree);
        assert(prob == mip);
        vector<int> left = *ref_left;
        vector<int> right = *ref_right;
        int cnt_vert = global_cnt_vert;
        int cnt_edj = ref_right->size();
        vector<vector<int> > edj(cnt_vert);
        for (int i = 0; i < cnt_edj; ++i) {
            if (glp_get_col_prim(prob, i + 1) > eps) {
                edj[left[i]].push_back(right[i]);
                edj[right[i]].push_back(left[i]);
            }
        }
        vector<char> visit(cnt_vert, 0);
        set<int> component;
        for (int i = 0; i < cnt_vert; ++i) {
            if (!visit[i]) {
                dfs(i, component, edj, visit);
                if (component.size() % 2 == 1) {
                    int half_size = component.size() / 2;
                    vector<double> ar(1, 0);
                    vector<int> ja(1, 0);
                    for (int k = 0; k < cnt_edj; ++k) {
                        if (component.find(left[k]) != component.end() &&
                            component.find(right[k]) != component.end()) {
                            ja.push_back(k + 1);
                            ar.push_back(1);

                        }
                    }
                    int idx = glp_add_rows(prob, 1);
                    glp_set_row_bnds(prob, idx, GLP_UP, half_size, half_size);
                    glp_set_mat_row(prob, idx, ja.size() - 1, ja.data(), ar.data());
                    return;
                }
                component.clear();
            }
        }
    }
    return;
}

int main() {
//    std::freopen("test/03", "r", stdin);
    cout << setprecision(6);
    cout << fixed;
    int cnt_vert, cnt_edj;
    cin >> cnt_vert >> cnt_edj;
    global_cnt_vert = cnt_vert;
    vector<int> left(cnt_edj);
    vector<int> right(cnt_edj);
    vector<int> weight(cnt_edj);
    for (int i = 0; i < cnt_edj; ++i) {
        cin >> left[i] >> right[i] >> weight[i];
    }
    ref_left = &left;
    ref_right = &right;
    glp_term_out(GLP_OFF);
    mip = glp_create_prob();
    glp_set_obj_dir(mip, GLP_MIN);
    glp_add_rows(mip, cnt_vert);
    for (int i = 0; i < cnt_vert; ++i) {
        glp_set_row_bnds(mip, i + 1, GLP_FX, 1, 1);
    }
    glp_add_cols(mip, cnt_edj);
    for (int i = 0; i < cnt_edj; ++i) {
        glp_set_col_kind(mip, i + 1, GLP_BV);
        glp_set_obj_coef(mip, i + 1, weight[i]);
    }
    vector<int> ia(1, 0);
    vector<int> ja(1, 0);
    vector<double> ar(1, 0);
    for (int i = 0; i < cnt_edj; ++i) {
        ia.push_back(left[i] + 1);
        ja.push_back(i + 1);
        ar.push_back(1);
        ia.push_back(right[i] + 1);
        ja.push_back(i + 1);
        ar.push_back(1);
    }
    glp_load_matrix(mip, static_cast<int>(ar.size()) - 1, ia.data(), ja.data(), ar.data());
    glp_smcp params_lp;
    glp_init_smcp(&params_lp);
    params_lp.msg_lev = GLP_MSG_OFF;
    glp_simplex(mip, &params_lp);
    glp_iocp param_mip;
    glp_init_iocp(&param_mip);
    param_mip.msg_lev = GLP_MSG_OFF;
    param_mip.br_tech = GLP_BR_MFV; /* most fractional variable */
    param_mip.bt_tech = GLP_BT_BLB; /* best local bound */
    param_mip.sr_heur = GLP_OFF; /* disable simple rounding heuristic */
    param_mip.gmi_cuts = GLP_ON; /* enable Gomory cuts */
    param_mip.cb_func = find_cut;
    glp_intopt(mip, &param_mip);
    int ans = glp_mip_obj_val(mip) + eps;
    cout << ans << '\n';
    for (int i = 0; i < cnt_edj; ++i) {
        if (glp_mip_col_val(mip, i + 1) > eps) {
            cout << i << ' ';
        }
    }
    cout << '\n';
    return 0;
}
