{ lib, pkgs, pythonPackages, inflection }:

pythonPackages.buildPythonPackage rec {
    pname = "diot";
    version = "0.1.4";
    name = "${pname}-${version}";

    src = pythonPackages.fetchPypi {
      inherit pname version;
      sha256 = "138583021f2b462f2d277d2c7dffb5eec8fa051f9a0337c71c1ead18aef7abb0";
    };

    propagatedBuildInputs = with pythonPackages; [
      inflection
    ];

    doCheck = false;

    meta = with lib; {
      homepage = https://github.com/pwwang/diot;
      description = "Python dictionary with dot notation.";
    };
}