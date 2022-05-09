{ lib, pkgs, pythonPackages, diot, pure-eval, varname, executing }:

pythonPackages.buildPythonPackage rec {
    pname = "pipda";
    version = "0.4.5";
    name = "${pname}-${version}";

    src = pythonPackages.fetchPypi {
      inherit pname version;
      sha256 = "0a2306427d2ce8c7a21fd6090d1148fe84c20623552ff5a76b1d8e94c6394d38";
    };

    propagatedBuildInputs = with pythonPackages; [
      diot
      pure-eval
      varname
      executing
    ];

    doCheck = false;

    meta = with lib; {
      homepage = "";
      description = "A framework for data piping in python";
    };
}
