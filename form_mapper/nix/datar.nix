{ lib, pkgs, pythonPackages, pipda }:


pythonPackages.buildPythonPackage rec {
    pname = "datar";
    version = "0.5.3";
    name = "${pname}-${version}";

    src = pythonPackages.fetchPypi {
      inherit pname version;
      sha256 = "b36f0be7570975d80ca0eec684822ed44e1f7030c514f6a89f5a4a5b592a56c7";
    };

    propagatedBuildInputs = with pythonPackages; [
      pipda
      pandas
    ];

    doCheck = false;

    meta = with lib; {
      homepage = https://github.com/pwwang/datar;
      description = "Port of dplyr and other related R packages in python, using pipda.";
    };
}
