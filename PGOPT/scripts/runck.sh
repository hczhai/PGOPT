#!/bin/bash

cat << EOF | $AFFCKHOME/bin/affck.exe > ck.out
{ parameter:
  ck_name = @GNAME;
  ck_number = @GNUMBER;
  ck_dimension = 3;
  ck_write_trajectory = false;
  ck_drel = @GDREL;
  ck_rotate_angle = 0;
  ck_output_xyz_files = pos_0.xyz;
  path_override = true;
  random_seed = @GRSEED;
}
{ exe: ck }
EOF
