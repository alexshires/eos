
eos-evaluate --kinematics s_min 1 --kinematics s_max 6 \
   --budget FF --vary "B->pi::f_+(0)@IKMvD-2014" \
    --observable "B->Pill::BR@LargeRecoil,l=mu"\
    | tee pimumu.low.value.data
   # --vary "B->pi::f_+(0)@BCL2008" --vary "B->pi::f_T(0)@BCL2008"\


