KPL/MK

   This meta-kernel lists a subset of kernels from the meta-kernel
   msgr_2011_v10.tm provided in the MESS-E/V/H-SPICE-6-V1.0 SPICE PDS3 archive,
   covering the whole or a part of the customer requested time period
   from 2011-08-12T21:40:06.000 to 2011-08-12T21:48:51.000.

   The documentation describing these kernels can be found in the
   complete MESS-E/V/H-SPICE-6-V1.0 SPICE PDS3 archive available at this URL

   https://naif.jpl.nasa.gov/pub/naif/pds/data/mess-e_v_h-spice-6-v1.0/messsp_1000

   To use this meta-kernel users may need to modify the value of the
   PATH_VALUES keyword to point to the actual location of the archive's
   ``data'' directory on their system. Replacing ``/'' with ``\''
   and converting line terminators to the format native to the user's
   system may also be required if this meta-kernel is to be used on a
   non-UNIX workstation.

   This meta-kernel was created by the NAIF node's SPICE PDS archive
   subsetting service version 2.1 on Wed Mar 29 02:57:05 PDT 2023.

 
   \begindata
 
      PATH_VALUES     = ( '/tmp/spice_kernels_messenger')
 
      PATH_SYMBOLS    = (
                         'KERNELS'
                        )
 
      KERNELS_TO_LOAD = (
                         '$KERNELS/lsk/naif0011.tls'
                         '$KERNELS/pck/pck00010_msgr_v23.tpc'
                         '$KERNELS/sclk/messenger_2548.tsc'
                         '$KERNELS/fk/msgr_dyn_v600.tf'
                         '$KERNELS/fk/msgr_v231.tf'
                         '$KERNELS/ik/msgr_epps_v100.ti'
                         '$KERNELS/ik/msgr_grns_v110.ti'
                         '$KERNELS/ik/msgr_mag_v021.ti'
                         '$KERNELS/ik/msgr_mascs_v100.ti'
                         '$KERNELS/ik/msgr_mdis_v160.ti'
                         '$KERNELS/ik/msgr_mla_v010.ti'
                         '$KERNELS/ik/msgr_rs_v111.ti'
                         '$KERNELS/ik/msgr_xrs_v001.ti'
                         '$KERNELS/spk/msgr_antenna_v000.bsp'
                         '$KERNELS/spk/msgr_040803_150430_150430_od431sc_2.bsp'
                         '$KERNELS/ck/msgr_1108_v02.bc'
                         '$KERNELS/ck/msgr_mdis_gm040819_150430v1.bc'

                         '$KERNELS/dsk/planets/mercury_m002_mes_v01.bds'
                        )
 
   \begintext
 

