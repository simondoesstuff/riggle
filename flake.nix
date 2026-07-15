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
    bits.url = "github:arq5x/bits";
    igd-src = {
      url = "github:databio/iGD";
      flake = false;
    };
    regioners-src = {
      url = "github:ACEnglish/regioners";
      flake = false;
    };
  };

  outputs =
    {
      nixpkgs,
      flake-utils,
      treefmt-nix,
      giggle,
      bits,
      igd-src,
      regioners-src,
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

        regioners = pkgs.rustPlatform.buildRustPackage {
          pname = "regioners";
          version = "0.3.1";
          src = regioners-src;
          cargoLock.lockFile = "${regioners-src}/Cargo.lock";
        };

        # for non-nix users https://hgdownload.gi.ucsc.edu/downloads.html#utilities_downloads
        liftOver = pkgs.stdenv.mkDerivation rec {
          pname = "liftOver";
          version = "latest";

          src =
            if pkgs.stdenv.isDarwin then
              pkgs.fetchurl {
                url = "https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/liftOver";
                sha256 = "sha256-e6uQT92vzWiK3WolABoXy0DiGqFtuwa4Z2W3Kc7kEOw=";
              }
            else
              pkgs.fetchurl {
                url = "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver";
                sha256 = "sha256-0000000000000000000000000000000000000000000="; # Placeholder
              };

          # Since we are fetching a bare binary, we need to skip the unpack phase
          dontUnpack = true;

          nativeBuildInputs = pkgs.lib.optionals pkgs.stdenv.isLinux [ pkgs.autoPatchelfHook ];
          buildInputs =
            pkgs.lib.optionals pkgs.stdenv.isLinux [
              pkgs.zlib
              pkgs.openssl
              pkgs.libpng
            ]
            ++ pkgs.lib.optionals pkgs.stdenv.isDarwin [ pkgs.xz ];

          installPhase =
            ''
              mkdir -p $out/bin
              cp $src $out/bin/liftOver
              chmod +x $out/bin/liftOver
            ''
            + pkgs.lib.optionalString pkgs.stdenv.isDarwin ''
              install_name_tool -change \
                /opt/homebrew/opt/xz/lib/liblzma.5.dylib \
                ${pkgs.xz.out}/lib/liblzma.5.dylib \
                $out/bin/liftOver
              codesign -s - $out/bin/liftOver 2>/dev/null || true
            '';
        };
      in
      {
        devShells.default = pkgs.mkShellNoCC {
          packages = with pkgs; [
            cargo
            # bio
            htslib # bgzip
            bedtools
            giggle.packages.${system}.default
            bits.packages.${system}.default
            igd
            liftOver
            regioners
            # misc
            just
            libiconv
            # analysis
            uv
            samply
            snakemake
          ];

          # Help the linker find libiconv on Darwin
          LIBRARY_PATH = pkgs.lib.optionalString pkgs.stdenv.isDarwin "${pkgs.libiconv}/lib";
        };

        formatter = treefmtEval.config.build.wrapper;
      }
    );
}
