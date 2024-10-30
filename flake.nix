{
    inputs = {
        nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
        scientific-fhs = {
            url = "github:manuelbb-upb/scientific-fhs/flake_module";
            inputs.nixpkgs.follows = "nixpkgs";
        };
    };
    
    outputs = inputs@{self, nixpkgs, scientific-fhs, ...}:
        let
            system = "x86_64-linux";
            pkgs = nixpkgs.legacyPackages.${system}.pkgs;
            julia-fhs = pkgs.callPackage scientific-fhs.fhsModule {
                enableJulia = true;
                enableConda = false;
                enableQuarto = false;
                enableNVIDIA = false;
                enableNode = false;
                juliaVersion = "1.10.4";
                commandScript = "julia";
                commandName = "julia";
            };
        in
        {
        devShells.${system}.default = pkgs.mkShell {
            buildInputs = [
                julia-fhs
                pkgs.gfortran      
            ];
            shellHook = ''
                export FREETYPE_ABSTRACTION_FONT_PATH="/run/current-system/sw/share/X11/fonts"
            '';
        };
    };
}
