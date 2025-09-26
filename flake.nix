{
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";

    scientific-nix-pkgs = {
      url = "github:manuelbb-upb/scientific-nix-pkgs";
      inputs.nixpkgs.follows = "nixpkgs";
    };

  };

  outputs = { 
    self, 
    nixpkgs, 
    scientific-nix-pkgs,
  }@inputs: 
  let
    system = "x86_64-linux"; 
    pkgs = nixpkgs.legacyPackages.${system}; 
    spkgs = scientific-nix-pkgs.packages.${system};

    julia = spkgs.julia-ld.override {
      version = "1.11.4";
      enable-matlab = false;
    };

  in
  {
    devShells.${system}.default = pkgs.mkShell {
      shellHook = ''
        export JULIA_NUM_THREADS=12
      '';
      packages = [
        julia
      ] ++ (with pkgs; [
        gfortran
      ]);
    };
  }; 
}
