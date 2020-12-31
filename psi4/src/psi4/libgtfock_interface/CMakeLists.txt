#list(APPEND sources
#  Minimal_Interface.cc
#  )
#
#psi4_add_module(lib gtfock_interface sources)

#target_include_directories(gtfock_interface
#  PUBLIC
#  )

include_directories([AFTER]
  /storage/home/hcoda1/2/jobrien60/p-cs207-0/testflight_data/psis/for_gtfock/psi4/psi4/src/psi4/libgtfock_interface/GTFock
  /storage/home/hcoda1/2/jobrien60/p-cs207-0/testflight_data/psis/for_gtfock/psi4/psi4/src/psi4/libgtfock_interface/GTFock/gtfock/pfock
  /storage/home/hcoda1/2/jobrien60/p-cs207-0/testflight_data/psis/for_gtfock/psi4/psi4/src/psi4/libgtfock_interface/GTFock/libcint
  /storage/home/hcoda1/2/jobrien60/p-cs207-0/testflight_data/psis/for_gtfock/psi4/psi4/src/psi4/libgtfock_interface/GTFock/GTMatrix
  /storage/home/hcoda1/2/jobrien60/.conda/envs/gtfock_JK/include
  )

