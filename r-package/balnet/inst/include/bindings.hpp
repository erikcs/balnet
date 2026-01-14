#pragma once
#include <Rcpp.h>
#include <RcppEigen.h>

#include <adelie_core/constraint/constraint_base.ipp>
#include <adelie_core/glm/glm_base.ipp>
#include <adelie_core/matrix/matrix_naive_base.ipp>
#include <adelie_core/matrix/matrix_naive_dense.ipp>
#include <adelie_core/state/state_gaussian_pin_base.ipp>
#include <adelie_core/state/state_gaussian_pin_naive.ipp>
#include <adelie_core/state/state_base.ipp>
#include <adelie_core/state/state_glm_naive.ipp>
#include <adelie_ext/glm_cbps.ipp>

namespace ad = adelie_core;

// Matrix type parameters
using value_t = double;
using index_t = int;
using bool_t = int;
using vec_value_t = ad::util::colvec_type<value_t>;
using vec_index_t = ad::util::colvec_type<index_t>;
using dense_64F_t = ad::util::colmat_type<value_t>;
using matrix_naive_dense_64F_t = ad::matrix::MatrixNaiveDense<ad::util::colmat_type<double>, int>;

// Constraint type parameters
using constraint_base_64_t = ad::constraint::ConstraintBase<double, int>;
using dyn_vec_constraint_t = std::vector<constraint_base_64_t*>;

// State type parameters
using vec_bool_t = ad::util::colvec_type<bool_t>;
using state_glm_naive_64_t = ad::state::StateGlmNaive<
    constraint_base_64_t, matrix_naive_dense_64F_t, double, int, int, int
>;

// Extension GLM
using glm_cbps_64_t = ad::glm::GlmCBPS<double>;


state_glm_naive_64_t* make_state(
    matrix_naive_dense_64F_t* X_ad,
    const Rcpp::List args
)
{
    const Eigen::Map<vec_value_t> eta = args["eta"];
    const Eigen::Map<vec_value_t> resid = args["resid"];
    const Rcpp::List constraints_r = args["constraints"];
    dyn_vec_constraint_t constraints;
    constraints.reserve(constraints_r.size());
    for (auto c : constraints_r) {
        constraints.push_back(nullptr);
    }
    const Eigen::Map<vec_index_t> groups = args["groups"];
    const Eigen::Map<vec_index_t> group_sizes = args["group_sizes"];
    const Eigen::Map<vec_index_t> dual_groups = args["dual_groups"];
    value_t alpha = args["alpha"];
    const Eigen::Map<vec_value_t> penalty = args["penalty"];
    const Eigen::Map<vec_value_t> offsets = args["offsets"];
    const Eigen::Map<vec_value_t> lmda_path = args["lmda_path"];
    value_t loss_null = args["loss_null"];
    value_t loss_full = args["loss_full"];
    value_t lmda_max = args["lmda_max"];
    value_t min_ratio = args["min_ratio"];
    size_t lmda_path_size = args["lmda_path_size"];
    size_t max_screen_size = args["max_screen_size"];
    size_t max_active_size = args["max_active_size"];
    value_t pivot_subset_ratio = args["pivot_subset_ratio"];
    size_t pivot_subset_min = args["pivot_subset_min"];
    value_t pivot_slack_ratio = args["pivot_slack_ratio"];
    const std::string screen_rule = args["screen_rule"];
    size_t irls_max_iters = args["irls_max_iters"];
    value_t irls_tol = args["irls_tol"];
    size_t max_iters = args["max_iters"];
    value_t tol = args["tol"];
    value_t adev_tol = args["adev_tol"];
    value_t ddev_tol = args["ddev_tol"];
    value_t newton_tol = args["newton_tol"];
    size_t newton_max_iters = args["newton_max_iters"];
    bool early_exit = args["early_exit"];
    bool setup_loss_null = args["setup_loss_null"];
    bool setup_lmda_max = args["setup_lmda_max"];
    bool setup_lmda_path = args["setup_lmda_path"];
    bool intercept = args["intercept"];
    size_t n_threads = args["n_threads"];
    const Eigen::Map<vec_index_t> screen_set = args["screen_set"];
    const Eigen::Map<vec_value_t> screen_beta = args["screen_beta"];
    const Eigen::Map<vec_bool_t> screen_is_active = args["screen_is_active"];
    size_t active_set_size = args["active_set_size"];
    const Eigen::Map<vec_index_t> active_set = args["active_set"];
    value_t beta0 = args["beta0"];
    value_t lmda = args["lmda"];
    const Eigen::Map<vec_value_t> grad = args["grad"];
    X_ad->mul(resid, vec_value_t::Ones(resid.size()), grad); // Init. grad to crossprod(X, resid)

    return new state_glm_naive_64_t(
        *X_ad, eta, resid, constraints, groups, group_sizes, dual_groups, alpha, penalty, offsets, lmda_path,
        loss_null, loss_full, lmda_max, min_ratio, lmda_path_size, max_screen_size, max_active_size,
        pivot_subset_ratio, pivot_subset_min, pivot_slack_ratio, screen_rule,
        irls_max_iters, irls_tol, max_iters, tol, adev_tol, ddev_tol, newton_tol, newton_max_iters,
        early_exit, setup_loss_null, setup_lmda_max, setup_lmda_path, intercept, n_threads,
        screen_set, screen_beta, screen_is_active, active_set_size, active_set, beta0, lmda, grad
    );
}

// Betas conversion utility, from adelie R's rcpp_state.cpp
template <class BetasType>
static auto convert_betas(
    size_t p,
    const BetasType& betas
)
{
    using value_t = typename std::decay_t<BetasType>::value_type::Scalar;
    using index_t = int;
    using vec_value_t = ad::util::rowvec_type<value_t>;
    using vec_index_t = ad::util::rowvec_type<index_t>;
    using sp_mat_t = Eigen::SparseMatrix<value_t, Eigen::RowMajor, index_t>;

    const size_t l = betas.size();
    size_t nnz = 0;
    for (const auto& beta : betas) {
        nnz += beta.nonZeros();
    }
    vec_value_t values(nnz);
    vec_index_t inners(nnz);
    vec_index_t outers(l+1);
    outers[0] = 0;
    int inner_idx = 0;
    for (size_t i = 0; i < l; ++i) {
        const auto& curr = betas[i];
        const auto nnz_curr = curr.nonZeros();
        Eigen::Map<vec_value_t>(
            values.data() + inner_idx,
            nnz_curr
        ) = Eigen::Map<const vec_value_t>(
            curr.valuePtr(),
            nnz_curr
        );
        Eigen::Map<vec_index_t>(
            inners.data() + inner_idx,
            nnz_curr
        ) = Eigen::Map<const vec_index_t>(
            curr.innerIndexPtr(),
            nnz_curr
        );
        outers[i+1] = outers[i] + nnz_curr;
        inner_idx += nnz_curr;
    }
    sp_mat_t out;
    out = Eigen::Map<const sp_mat_t>(
        l,
        p,
        nnz,
        outers.data(),
        inners.data(),
        values.data()
    );
    return out;
}
