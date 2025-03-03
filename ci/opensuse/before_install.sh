zypper update -y \
    && zypper install -y \
              --oldpackage glibc=2.40 glibc-locale-base=2.40 \
              cmake cgal-devel ninja gcc gcc-c++ \
              libboost_serialization-devel \
              libboost_timer-devel \
              libboost_program_options-devel \
              libboost_filesystem-devel \
              libboost_test-devel \
    && rm -rf build
