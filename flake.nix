{
  description = "Rust Giggle";
  inputs = {
    # nix stuff
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.11";
    flake-utils.url = "github:numtide/flake-utils";
    treefmt-nix = {
      url = "github:numtide/treefmt-nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    # packages
    giggle.url = "path:/Users/simon/Code/lab/giggle-dev/giggle";
    igd-src = {
      url = "github:databio/iGD";
      flake = false;
    };
  };

  outputs =
    {
      nixpkgs,
      flake-utils,
      treefmt-nix,
      giggle,
      igd-src,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};

        treefmtEval = treefmt-nix.lib.evalModule pkgs {
          projectRootFile = "flake.nix";

          # Enable the formatters you need
          programs.nixpkgs-fmt.enable = true; # Formats .nix files
          programs.rustfmt.enable = true; # Formats .rs files via cargo
        };

        igd = pkgs.stdenv.mkDerivation {
          pname = "igd";
          version = "unstable";
          src = igd-src;
          buildInputs = [ pkgs.zlib ];
          buildPhase = "make";
          installPhase = ''
            mkdir -p $out/bin
            cp bin/igd $out/bin/
          '';
        };
      in
      {
        devShells.default = pkgs.mkShellNoCC {
          packages = with pkgs; [
            just
            htslib # bgzip
            bedtools
            giggle.packages.${system}.default
            libiconv
            cargo
            samply
            igd
          ];

          # Help the linker find libiconv on Darwin
          LIBRARY_PATH = pkgs.lib.optionalString pkgs.stdenv.isDarwin "${pkgs.libiconv}/lib";
        };

        formatter = treefmtEval.config.build.wrapper;
      }
    );
}
