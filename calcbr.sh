#eos-evaluate --range s 14.18 22.86 12\
#    --observable "B->K^*ll::BR@LargeRecoil,form-factors=KMPW2010,l=mu"\
    #--kinematics s_min 1 --kinematics s_max 6 \
    #--budget CKM --vary CKM::lambda --vary CKM::A --vary CKM::rhobar --vary CKM::etabar 
    #--budget FF --vary "B->K^*::F^V(0)@KMPW2010"\


#eos-evaluate --range s_min 14.18 s_max 22.86 12\
#    --observable "B->K^*ll::BR@LargeRecoil,form-factors=KMPW2010,l=e"\
#--kinematics s_min 1 --kinematics s_max 6 \
#--budget "FF" --vary "B->K^*::F^V(0)@KMPW2010"

eos-evaluate --kinematics s_min 0.1 --kinematics s_max 10 \
    --observable "B->Pill::BR@LargeRecoil,l=mu"\

eos-evaluate --kinematics s_min 10 --kinematics s_max 26 \
    --observable "B->Pill::BR@LowRecoil,l=mu"\

eos-evaluate --kinematics s_min 0.1 --kinematics s_max 10 \
    --observable "B->Kll::BR@LargeRecoil,l=mu"\

eos-evaluate --kinematics s_min 14 --kinematics s_max 22 \
    --observable "B->Kll::BR@LowRecoil,l=mu"\






eos-evaluate --kinematics s_min 1 --kinematics s_max 6 \
   --budget FF  \
    --vary "B->pi::f_+(0)@BCL2008" --vary "B->pi::f_T(0)@BCL2008"\
    --observable "B->Pill::BR@LargeRecoil,l=mu"\
    | tee pimumu.low.value.data

eos-evaluate --kinematics s_min 15 --kinematics s_max 20 \
   --budget FF --vary "B->pi::f_+(0)@BCL2008" --vary "B->pi::f_T(0)@BCL2008"\
    --observable "B->Pill::BR@LowRecoil,l=mu"\
    | tee pimumu.large.value.data

eos-evaluate --range s 0.1 9 60 \
   --budget FF --vary "B->pi::f_+(0)@BCL2008" --vary "B->pi::f_T(0)@BCL2008"\
    --observable "B->Pill::dBR/ds@LargeRecoil,l=mu"\
    > pimumu.large.data

eos-evaluate --range s 12 25 50 \
    --observable "B->Pill::dBR/ds@LowRecoil,l=mu" \
   --budget FF --vary "B->pi::f_+(0)@BCL2008" --vary "B->pi::f_T(0)@BCL2008"\
    > pimumu.low.data
#
eos-evaluate --kinematics s_min 1 --kinematics s_max 6 \
    --budget FF --vary "B->K::F^p(0)@KMPW2010" --vary "B->K::F^t(0)@KMPW2010"\
    --observable "B->Kll::BR@LargeRecoil,l=mu"\
    | tee kmumu.low.value.data

eos-evaluate --kinematics s_min 15 --kinematics s_max 20 \
    --budget FF --vary "B->K::F^p(0)@KMPW2010" --vary "B->K::F^t(0)@KMPW2010"\
    --observable "B->Kll::BR@LowRecoil,l=mu"\
    | tee kmumu.large.value.data

eos-evaluate --range s 0.1 9 60 \
    --observable "B->Kll::dBR/ds@LargeRecoil,l=mu"\
    --budget FF --vary "B->K::F^p(0)@KMPW2010" --vary "B->K::F^t(0)@KMPW2010"\
    > kmumu.large.data

eos-evaluate --range s 12 25 50 \
    --observable "B->Kll::dBR/ds@LowRecoil,l=mu" \
    --budget FF --vary "B->K::F^p(0)@KMPW2010" --vary "B->K::F^t(0)@KMPW2010"\
    > kmumu.low.data

    
