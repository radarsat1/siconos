
// List of solvers generated with
//   git grep '^ *SICONOS_[A-Z_]* = [0-9]' | sed 's/.*\(SICONOS_[A-Z_]*\).*/SICONOS_SOLVER_MACRO(\1); \\/' |
//   grep -v NUMERICS_PROBLEM > NonSmoothSolvers/SiconosNumerics_Solvers.h


#undef SICONOS_SOLVER_MACRO
#define SICONOS_REGISTER_SOLVERS() \
SICONOS_SOLVER_MACRO(SICONOS_AVI_CAOFERRIS); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_LEMKE); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_NSGS_SBM); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_PGS); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_CPG); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_LATIN); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_LATIN_W); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_QP); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_NSQP); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_NEWTONMIN); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_NEWTON_FBLSA); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_PSOR); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_RPGS); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_PATH); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_ENUM); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_AVI_CAOFERRIS); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_PIVOT); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_BARD); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_MURTY); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_NEWTON_MINFBLSA); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_PATHSEARCH); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_PIVOT_LUMOD); \
SICONOS_SOLVER_MACRO(SICONOS_LCP_GAMS); \
SICONOS_SOLVER_MACRO(SICONOS_MCP_FB); \
SICONOS_SOLVER_MACRO(SICONOS_MCP_NEWTON_FBLSA); \
SICONOS_SOLVER_MACRO(SICONOS_MCP_NEWTON_MINFBLSA); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_PGS); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_RPGS); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_PSOR); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_RPSOR); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_PATH); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_ENUM); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_SIMPLEX); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_DIRECT_ENUM); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_PATH_ENUM); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_DIRECT_SIMPLEX); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_DIRECT_PATH); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_DIRECT_PATH_ENUM); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_FB); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_DIRECT_FB); \
SICONOS_SOLVER_MACRO(SICONOS_MLCP_PGS_SBM); \
SICONOS_SOLVER_MACRO(SICONOS_NCP_NEWTON_FBLSA); \
SICONOS_SOLVER_MACRO(SICONOS_NCP_NEWTON_MINFBLSA); \
SICONOS_SOLVER_MACRO(SICONOS_NCP_PATHSEARCH); \
SICONOS_SOLVER_MACRO(SICONOS_NCP_PATH); \
SICONOS_SOLVER_MACRO(SICONOS_RELAY_PGS); \
SICONOS_SOLVER_MACRO(SICONOS_RELAY_ENUM); \
SICONOS_SOLVER_MACRO(SICONOS_RELAY_PATH); \
SICONOS_SOLVER_MACRO(SICONOS_RELAY_LEMKE); \
SICONOS_SOLVER_MACRO(SICONOS_RELAY_LATIN); \
SICONOS_SOLVER_MACRO(SICONOS_RELAY_NLGS); \
SICONOS_SOLVER_MACRO(SICONOS_RELAY_AVI_CAOFERRIS); \
SICONOS_SOLVER_MACRO(SICONOS_RELAY_AVI_CAOFERRIS_TEST); \
SICONOS_SOLVER_MACRO(SICONOS_VI_EG); \
SICONOS_SOLVER_MACRO(SICONOS_VI_FPP); \
SICONOS_SOLVER_MACRO(SICONOS_VI_HP); \
SICONOS_SOLVER_MACRO(SICONOS_VI_BOX_QI); \
SICONOS_SOLVER_MACRO(SICONOS_VI_BOX_AVI_LSA);
