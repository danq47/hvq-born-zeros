c -*- Fortran -*-

c flg_nlotest: perform calls to outfun for testing NLO output
c flg_withsubtr: counterterms are included in NLO test
c flg_withdamp: Born zero procedure
c flg_withreg: if regular regions exist or not
c flg_smartsig: remember or not equal squared amplitude
c flg_in_smartsig: set to true when the program is comuputing
c                  amplitudes to determine identical ones.
c                  Can be use to turn off kinematical cuts
c                  in the amplitudes
c flg_bornonly: do the Born contribution only
c flg_LOevents: generate events using only the Born. Real and virtual contributions
c               are identically zero

c flg_lightpart_check: in genflavreglist perform or not perform the check
c                     that there are no coloured light partons before flst_lightpart

c flg_debug: outputs extra infos for radiation region and randon numbers on the LHEF   
c flg_withnegweights: allows negative weighted events
c flg_jacsing:  importance sampling for singular jacobians in case of massless FSR
c flg_weightedev:  allows weighted events
c flg_pdfreweight: outputs extra infos useful for pdf's reweighting on the LHEF
c flg_collremnsamp:  importance sampling for collinear remnants
c flg_reweight: outputs extra infos for reweighting LH events

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C DQ flags added
C
C flg_newsuda: .true. if using the new Sudakov to generate
C radiation, and .false. if using the old Sudakov
C 
C flg_remnant: .true. if the event is a remnant,
C and .false. otherwise
C
C flg_b0_NLO: .true. - testing Born zeros.
C
C flg_rsoft: .true. means we use the soft limit of the planar
C matrix elements, .false. means we use the usual (hard) MEs
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      logical flg_nlotest,flg_withsubtr,flg_withdamp,flg_withreg,
     1     flg_smartsig,flg_in_smartsig,flg_bornonly,flg_debug,
     2     flg_withnegweights,
     3     flg_jacsing,flg_weightedev,flg_pdfreweight,flg_collremnsamp,
     4     flg_lightpart_check,flg_btlscalereal,flg_btlscalect,
     5     flg_bornzerodamp,flg_ckkwscalup,
     6     flg_minlo,flg_minlo_nnll,flg_minlo_real,
     7     flg_reweight,flg_newweight,flg_fastbtlbound,
     8     flg_storemintupb,flg_doublefsr,flg_monitorubound,
     9     flg_drscheme,flg_withresrad,flg_with_em,flg_em_rad,
     $     flg_LOevents,flg_evenmaxrat,flg_novirtual,flg_noevents,
     1     flg_doubletags,flg_analysisextrainfo,flg_fullrwgt,
C
C DQ flags
C
     2     flg_newsuda,flg_remnant,flg_b0_NLO,flg_rsoft,
C
C End DQ flags
C
     3     flg_ubexcess_correct,flg_mintupb_xless
      character * 1 flg_btildepart
      character * 20 flg_processid
      common/pwhg_flg/flg_nlotest,flg_withsubtr,flg_withdamp,
     1     flg_withreg,flg_smartsig,flg_in_smartsig,
     2     flg_bornonly,flg_debug,flg_LOevents,
     3     flg_withnegweights,flg_jacsing,flg_weightedev,
     4     flg_pdfreweight,flg_collremnsamp,flg_lightpart_check,
     5     flg_btlscalereal,flg_btlscalect,
     6     flg_bornzerodamp,flg_ckkwscalup,
     7     flg_minlo,flg_minlo_nnll,flg_minlo_real,
     8     flg_reweight,flg_newweight,flg_fastbtlbound,
     9     flg_storemintupb,flg_doublefsr,flg_monitorubound,
     1     flg_drscheme,flg_withresrad,flg_with_em,flg_em_rad,
     2     flg_evenmaxrat,flg_novirtual,flg_noevents,flg_doubletags,
     3     flg_analysisextrainfo,flg_fullrwgt,flg_ubexcess_correct,
C
C DQ flags
C
     4     flg_newsuda,flg_remnant,flg_b0_NLO,flg_rsoft,
C
C End DQ flags
C
     5     flg_mintupb_xless,
c     cstrings      
     2     flg_btildepart,flg_processid
      save /pwhg_flg/
