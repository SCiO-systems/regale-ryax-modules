[(
  self: super:
  let
      myPythonPackages = super.callPackage ./nix { };
      packageOverrides = python-self: python-super: myPythonPackages;
  in {
      python3 = super.python3.override {inherit packageOverrides;};
    }
)]
