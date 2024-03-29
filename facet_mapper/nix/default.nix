{ lib, pkgs, python, python3Packages }:

let
  self = rec {
    callPackage = lib.callPackageWith (pkgs // python3Packages // self);

    # The folowing comment line is use to inject packages, do not change it.
    # ADD NEW PACKAGE HERE
    # executing = callPackage ./executing.nix { };
    varname = callPackage ./varname.nix { };
    diot = callPackage ./diot.nix { };
    pipda = callPackage ./pipda.nix { };
    datar = callPackage ./datar.nix { };
  };
in self
