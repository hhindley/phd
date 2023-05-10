
module ModelsModule
export rtc_model, rtc_all_t!, atp_gr_t!, atp_kin_t!, lam_kin_t, atp_t!, rtc_model_kin, rtc_model_lam_t!
include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl")
include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl")
include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl")
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl")
end