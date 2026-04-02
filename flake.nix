{
  description = "Rust Giggle";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.05";
    flake-utils.url = "github:numtide/flake-utils";
    giggle.url = "path:/Users/simon/Code/lab/giggle-dev/giggle";
  };

  outputs =
    {
      nixpkgs,
      flake-utils,
      giggle,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
      {
        devShells.default = pkgs.mkShellNoCC {
          packages = with pkgs; [
            just
            htslib # bgzip
            giggle.packages.${system}.default
            cargo
          ];
        };
      }
    );
}
