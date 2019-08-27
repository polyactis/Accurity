extern crate cc;

/*
use std::env;
use std::fs;
use std::path::PathBuf;
*/


fn main() {
    //let out = PathBuf::from(env::var_os("OUT_DIR").unwrap());
    //fs::remove_dir_all(&out).unwrap();
    //fs::create_dir(&out).unwrap();

    cc::Build::new()
        .cpp(true);
        /*
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-lgsl")
        .flag_if_supported("-lgslcblas")
        .file("Accurity/infer_class.cpp")
        .file("Accurity/BaseGADA.cc")
        .file("Accurity/read_para.cpp")
        .file("Accurity/prob.cpp")
        .compile("libaccurity.a");
        */
}