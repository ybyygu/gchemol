
# Table of Contents

1.  [Features](#orgebbb3d6)
2.  [Installing](#org374e644)
    1.  [install rust:](#org9f8c853)
    2.  [add Cargo.toml in your project:](#orge65a875)
3.  [Molecule](#orgf78138d)
4.  [Read/write chemical files](#org02078bd)
5.  [Templating](#orgb084e82)
6.  [Related projects](#orgc6b83b7)

gchemol is a graph-based chemical object library implemented in Rust programming
language.

[![Build Status](https://travis-ci.org/ybyygu/gchemol.svg?branch=master)](https://travis-ci.org/ybyygu/gchemol)
[![GPL3 licensed](https://img.shields.io/badge/license-GPL3-blue.svg)](./LICENSE)


<a id="orgebbb3d6"></a>

# Features

-   Fast and safe.
-   Easy to deploy in server environment.
-   core graph data structure using [petgraph](https://github.com/bluss/petgraph)
-   read molecules in various formats using [nom](https://github.com/Geal/nom) parser combinators
-   linear algebra backed by [nalgebra](http://nalgebra.org/)
-   render molecule in user defined formats by templating with [handlebars](https://github.com/sunng87/handlebars-rust)


<a id="org374e644"></a>

# Installing


<a id="org9f8c853"></a>

## install rust:

TBD


<a id="orge65a875"></a>

## add Cargo.toml in your project:

    [dependencies]
    gchemol = {git = "https://github.com/ybyygu/gchemol"}


<a id="orgf78138d"></a>

# Molecule

TBD


<a id="org02078bd"></a>

# Read/write chemical files

TBD


<a id="orgb084e82"></a>

# Templating

For syntax details: [sunng87/handlebars-rust: Rust templating with Handlebars](https://github.com/sunng87/handlebars-rust)

molecule title

    {{molecule_title}}

total number of atoms in molecule

    {{molecule.number_of_atoms}}

total number of bonds in molecule

    {{molecule.number_of_bonds}}

loop over atoms

    {{#each molecule.atoms as |a| ~}}
    {{a.symbol}} {{a.x}} {{a.y}} {{a.z}}
    {{/each~}}

atom properties

    {{atom.symbol}}

cartesian coordinates

    {{atom.x}} {{atom.y}} {{atom.z}}

fractional coordinates

    {{atom.fx}} {{atom.fy}} {{atom.fz}}

unit cell parameters

    {{molecule.unit_cell.a}}
    {{molecule.unit_cell.b}}
    {{molecule.unit_cell.c}}
    {{molecule.unit_cell.alpha}}
    {{molecule.unit_cell.beta}}
    {{molecule.unit_cell.gamma}}
    {{molecule.unit_cell.va}}
    {{molecule.unit_cell.vb}}
    {{molecule.unit_cell.vc}}

format string or number:

    {{format 12.16 width=12 prec=6 align="right"}}

alignment control:

    {{align="left"}}
    {{align="right"}}
    {{align="center"}}


<a id="orgc6b83b7"></a>

# Related projects

-   lumol
-   ase
-   pymatgen

