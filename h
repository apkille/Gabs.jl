[33mcommit d480e0b22cdbe6d3518c1219a4f4b52ca23251d1[m[33m ([m[1;36mHEAD -> [m[1;32mfa/i39[m[33m)[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Tue Mar 11 11:00:36 2025 +0500

    enhance tests in test_randoms.jl

[33mcommit 7cb4314120e6dcb3bc81b4177a2ae0bedbf77fad[m[33m ([m[1;31morigin/fa/i39[m[33m)[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Tue Mar 11 10:49:45 2025 +0500

    simplify the type-based dispatch

[33mcommit b5e53c506f8d24750248e8ba709e94b60918b617[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Tue Mar 11 10:26:43 2025 +0500

    fixup

[33mcommit 8f660b6bf05e9dbe3812022dff89736edede82c5[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Tue Mar 11 10:11:54 2025 +0500

    improvements

[33mcommit 631a001e45ab8a606928c04d8a7c75809264f825[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Tue Mar 11 10:10:05 2025 +0500

    improvements

[33mcommit 852647915a126dfe76a2dba5132e8117a0723f6c[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Tue Mar 11 10:02:12 2025 +0500

    improvements

[33mcommit 1ced61e063ddfac45e9ffd1113136660174c886e[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Sun Mar 9 15:26:14 2025 +0500

    cleanup

[33mcommit 9c27458ba891e335c359df0b4ceb2914e9e1ae23[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Sun Mar 9 15:12:16 2025 +0500

    improvements

[33mcommit b4dc7c84b0fa503a3ab506c0f2a6fae10acf8b86[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 17:37:31 2025 +0500

    undo tests that used previous notation of using static arrays

[33mcommit e10368bb625e6fbfcc95816932aba777d6fbf17e[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 17:32:37 2025 +0500

    polish tests

[33mcommit 920f21673b8c379c90fc161f8b658e5eb51c389d[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 16:05:52 2025 +0500

    polish

[33mcommit 4fee0c8dd9b44c5664b135caba0a71ce2cd0ea32[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 15:46:34 2025 +0500

    polish

[33mcommit 68e5665044d79c9408a47cc063f65ef8e80b53e4[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 15:27:47 2025 +0500

    polish

[33mcommit ea6dc66e4902628f16d9c7dd9f3d80452cf6c49f[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 15:25:22 2025 +0500

    add wonderful codereview suggestions: improvements and polish

[33mcommit fd7050e829396c5a1b901748cfaff4d26f9954f3[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 15:13:01 2025 +0500

    use dtype and ttype :)

[33mcommit 2594d9739b08c7dd05ec6bc0e76b786a3a573651[m
Author: Feroz Ahmad <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 15:10:14 2025 +0500

    Update src/channels.jl
    
    Co-authored-by: Andrew Kille <68079167+apkille@users.noreply.github.com>

[33mcommit f8ec8bc6a3a83e3aa9a8956172dc9faaf23f870f[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 15:09:46 2025 +0500

    add test broken error that seems to be due to equality checks for Num type in Symbolics that we previously encountered

[33mcommit 67f594bd510e88f8152a6a0f871f3450ba1f05e5[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 15:02:39 2025 +0500

    add test_broken for some errors and comment out tests for phaseshift as they throw error

[33mcommit 07e20fc13a170924e26ea67dc0debe54ffb37f66[m
Merge: 589be9c e1c8a19
Author: Feroz Ahmad <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 14:54:10 2025 +0500

    Merge branch 'apkille:main' into fa/i39

[33mcommit 589be9cf2a9c203e040a1cdbde4d7027bece8c83[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 14:43:41 2025 +0500

    add wonderful codereview suggestions: improvements to type-based dispatch

[33mcommit ec9cd51ca4f264afcc3d963050d5d2de872e4cda[m
Author: Feroz Ahmad <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 14:39:00 2025 +0500

    Update ext/StaticArraysExt/utils.jl
    
    Co-authored-by: Andrew Kille <68079167+apkille@users.noreply.github.com>

[33mcommit e1c8a19704187a292b2fa9c88f1f840465c37e4e[m
Author: Feroz Ahmad <feroz.ahmad.email@gmail.com>
Date:   Thu Feb 27 14:38:16 2025 +0500

    validate symbolic Gaussian channels (#50)
    
    * validate symbolic channels
    
    * add Symbolics to test suite
    
    * add symbolic displacement and improve tests
    
    * improve tests
    
    * polish
    
    * polish
    
    * add wonderful codereview suggestions
    
    * add wonderful codereview suggestion: simplify tests
    
    * add wonderful codereview suggestion: simplify tests
    
    * add wonderful codereview suggestions: simplify tests
    
    * add wonderful codereview suggestions: using isapprox as == throws error
    
    * add wonderful codereview suggestions: simplify tests
    
    * add codereview wonderful suggestions: improvements
    
    * polish
    
    * Add ħ = 2 as default Gabs convention and clean up docstrings (#49)
    
    * add ħ = 2 as default Gabs convention
    
    * add error message
    
    * add more hbar checks
    
    * fix docstring errors
    
    * fix docstrings and add docs about hbar convention
    
    * add more hbar tests
    
    * Validate symbolic Gaussian unitaries (#46)
    
    * Validate Symbolic Unitaries
    
    * fixup
    
    * add Symbolic two-mode squeeze tests
    
    * Symbolic phase-shift operator tests
    
    * combine symbolic beamsplitter and phase-shift operator tests into one test
    
    * combine single and two-mode squeeze to one testset to avoid redundancy
    
    * Symbolic tensor products for gaussian unitaries -  tests
    
    * polish tests
    
    * some tweaks to enable symbolic action and add basic tests
    
    * undo tweaks
    
    * add codereview suggestions
    
    * add codereview suggestion: basic symbolic action tests
    
    * add test_broken for symbolic action tests
    
    * improve tests by adding apply on squeezedstate
    
    * polish
    
    * fix merge conflicts
    
    * fix squeezed state sign
    
    * update changelog
    
    ---------
    
    Co-authored-by: Feroz Ahmad <feroz.ahmad.email@gmail.com>
    
    * cleanup tests
    
    * fix indexing issue by ensuring r and n are vectors in _amplifier
    
    * add wonderful suggestions: symbolic tensor product and symbolic actions for channels
    
    * define alpha2
    
    * clean up tests
    
    * address amplifier typo
    
    ---------
    
    Co-authored-by: Andrew Kille <68079167+apkille@users.noreply.github.com>
    Co-authored-by: apkille <killeandrew@gmail.com>

[33mcommit c6f2f83f1bb65dac38de441db845b1c15506094f[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Wed Feb 26 15:02:30 2025 +0500

    add missing comments that got removed

[33mcommit f33ef1404b85e5f7adc256eefe4a714ddd881594[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Wed Feb 26 14:51:56 2025 +0500

    add more tests for SArray as well

[33mcommit 783ec5ef911db72886fd7675e87dfefba40bd4ec[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Wed Feb 26 14:38:09 2025 +0500

    add wonderful codereview suggestions: polish

[33mcommit df7612e1a9eec3b98691afd29c01656106858249[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Wed Feb 26 13:48:05 2025 +0500

    add wonderful codereview suggestions: polish

[33mcommit 9db705d4c48315ad439bca28ba8668076bfa2dc9[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Wed Feb 26 13:38:33 2025 +0500

    add wonderful codereview suggestions: minor fixes/generalizations to _tensor and _promote_output_matrix

[33mcommit a7ba8223e2bd52c62f93fce5a8d5ab15cdbf5fdf[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Wed Feb 26 12:29:02 2025 +0500

    add codereview suggestions: use _infer_types via use traits instead of dispatch

[33mcommit f85baa33183581fba3c31ea59356096b4e2ed60e[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Sat Feb 22 17:24:15 2025 +0500

    fixup to resolve merge conflicts errors

[33mcommit dea74cd31af8b42375c0125c95cb3abeaf56b339[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Sat Feb 22 16:51:45 2025 +0500

    fix errors after merge conflict

[33mcommit 6721955b869815f0b302f7c514ffbb731c229101[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Sat Feb 22 16:36:37 2025 +0500

    fix merge conflict typo caused by resolving merge conflict

[33mcommit cd8c8e2e84623a51cffe6ead35feecbdf35e1857[m
Merge: 1172889 d6b1330
Author: Feroz Ahmad <feroz.ahmad.email@gmail.com>
Date:   Sat Feb 22 16:29:27 2025 +0500

    Merge branch 'main' into fa/i39

[33mcommit d6b133058ffea511dad09fffc3babe380c9a30d5[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Fri Feb 21 23:46:36 2025 -0500

    Add ħ = 2 as default Gabs convention and clean up docstrings (#49)
    
    * add ħ = 2 as default Gabs convention
    
    * add error message
    
    * add more hbar checks
    
    * fix docstring errors
    
    * fix docstrings and add docs about hbar convention
    
    * add more hbar tests
    
    * Validate symbolic Gaussian unitaries (#46)
    
    * Validate Symbolic Unitaries
    
    * fixup
    
    * add Symbolic two-mode squeeze tests
    
    * Symbolic phase-shift operator tests
    
    * combine symbolic beamsplitter and phase-shift operator tests into one test
    
    * combine single and two-mode squeeze to one testset to avoid redundancy
    
    * Symbolic tensor products for gaussian unitaries -  tests
    
    * polish tests
    
    * some tweaks to enable symbolic action and add basic tests
    
    * undo tweaks
    
    * add codereview suggestions
    
    * add codereview suggestion: basic symbolic action tests
    
    * add test_broken for symbolic action tests
    
    * improve tests by adding apply on squeezedstate
    
    * polish
    
    * fix merge conflicts
    
    * fix squeezed state sign
    
    * update changelog
    
    ---------
    
    Co-authored-by: Feroz Ahmad <feroz.ahmad.email@gmail.com>

[33mcommit 1172889f7513641fdb7970a5e282b749ed999961[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Fri Feb 21 22:33:32 2025 +0500

    cleaner dispatching for gaussian channels

[33mcommit ca20715e5d92e88b992c394b4deb7f4e41ca76e0[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Fri Feb 21 22:03:45 2025 +0500

    use cleaner dispatching onto static arrays for gaussian unitaries

[33mcommit 5b4246050f03033aa898caf7e95b8d1dd15dd8ef[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Fri Feb 21 18:48:55 2025 +0500

    some tests are broken

[33mcommit 964cb43e28cd23da271ab4df2f485b968f82a485[m
Author: Fe-r-oz <feroz.ahmad.email@gmail.com>
Date:   Fri Feb 21 18:40:27 2025 +0500

     fix #39 Cleaner dispatching onto StaticArrays

[33mcommit 7204e9e62825d7562bd7137d6d148ac990a46504[m
Author: Feroz Ahmad <feroz.ahmad.email@gmail.com>
Date:   Fri Feb 21 05:07:47 2025 +0500

    Validate symbolic Gaussian unitaries (#46)
    
    * Validate Symbolic Unitaries
    
    * fixup
    
    * add Symbolic two-mode squeeze tests
    
    * Symbolic phase-shift operator tests
    
    * combine symbolic beamsplitter and phase-shift operator tests into one test
    
    * combine single and two-mode squeeze to one testset to avoid redundancy
    
    * Symbolic tensor products for gaussian unitaries -  tests
    
    * polish tests
    
    * some tweaks to enable symbolic action and add basic tests
    
    * undo tweaks
    
    * add codereview suggestions
    
    * add codereview suggestion: basic symbolic action tests
    
    * add test_broken for symbolic action tests
    
    * improve tests by adding apply on squeezedstate
    
    * polish

[33mcommit bf5fdc3db07195d9243e7ecdb7dc03b3b6839ca0[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Thu Feb 13 22:23:55 2025 -0500

    Add seminal papers section to documentation (#48)
    
    * add seminal papers to docs
    
    * update changelog and project.toml
    
    * extend "ome" typo

[33mcommit 541ae0b4fb94c58e009268ecfa44d5c7a727370b[m
Author: Feroz Ahmad <feroz.ahmad.email@gmail.com>
Date:   Tue Feb 11 09:13:41 2025 +0500

    fix: correct the rendering of Latexify-based covariance matrix (#47)
    
    * fix: correct the rendering of Latexify-based covariance matrix
    
    * Update docs/src/tutorials.md
    
    Co-authored-by: Andrew Kille <68079167+apkille@users.noreply.github.com>
    
    ---------
    
    Co-authored-by: Andrew Kille <68079167+apkille@users.noreply.github.com>

[33mcommit 332c27e94691d126a23e3e83f1a07cacb8cfd40c[m
Author: apkille <killeandrew@gmail.com>
Date:   Thu Feb 6 20:50:36 2025 -0500

    update changelog and bump to v1.2.8

[33mcommit 4c4956e3b89248f34010a93e57160d2f90bf0858[m
Author: Feroz Ahmad <feroz.ahmad.email@gmail.com>
Date:   Fri Feb 7 06:25:01 2025 +0500

    fix #43: Validate and add docs for Gaussian states containing symbolic variables (#44)
    
    * resolve #43 Validate and add docs for Gaussian objects containing symbolic variables
    
    * fix extra space
    
    * Symbolic squeezed states
    
    * Symbolic coherent states
    
    * remove stupid mistake - add Symbolics to tests and docs, Latexify to docs
    
    * add Missings to typos
    
    * remove Latexify as dependency - just use it in docs via doc project toml
    
    * polish
    
    * shift documentation from manual to tutorial
    
    * fix small typo
    
    * Update docs/src/tutorials.md
    
    Co-authored-by: Andrew Kille <68079167+apkille@users.noreply.github.com>
    
    * Update docs/src/tutorials.md
    
    Co-authored-by: Andrew Kille <68079167+apkille@users.noreply.github.com>
    
    * Update docs/src/tutorials.md
    
    Co-authored-by: Andrew Kille <68079167+apkille@users.noreply.github.com>
    
    * add codereview suggestions: improve tutorial
    
    * add codereview suggestions: use rs and thetas for squeezedstates
    
    * add code review suggestion: alphas for coherent state
    
    * add code review suggestion: remove redundancy in EPR tests
    
    * symbolic tensor products :)
    
    * fix the doctest error that was causes by using Latexify in doctest
    
    * Symbolic thermal states and tests
    
    * undo a change to project toml
    
    * use float(eltype(P) so thermal states work with symbolic and ints etc
    
    * Update docs/src/tutorials.md
    
    Co-authored-by: Andrew Kille <68079167+apkille@users.noreply.github.com>
    
    * Update docs/src/tutorials.md
    
    Co-authored-by: Andrew Kille <68079167+apkille@users.noreply.github.com>
    
    * Update docs/make.jl
    
    Co-authored-by: Andrew Kille <68079167+apkille@users.noreply.github.com>
    
    * Create .gitignore
    
    Add manifest in .gitignore
    
    * rm the first test set
    
    * add some tests from removed testset 1 to testset2
    
    * don't show output of newst.covar in doctest, we are already doing it afterwards
    
    * polish testset 1
    
    * Delete Manifest.toml
    
    ---------
    
    Co-authored-by: Andrew Kille <68079167+apkille@users.noreply.github.com>

[33mcommit e5d5398b30e09262535c7a2571c8c6d452ec221c[m[33m ([m[1;31morigin/main[m[33m, [m[1;31morigin/HEAD[m[33m, [m[1;32mmain[m[33m)[m
Author: apkille <killeandrew@gmail.com>
Date:   Thu Jan 30 19:01:23 2025 -0500

    bump SympFact to v0.1.5

[33mcommit 303e71aabec5d2667707ea0478fb3c7e2f8b57b9[m
Author: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>
Date:   Mon Jan 27 19:30:32 2025 -0500

    Bump dawidd6/action-download-artifact from 7 to 8 (#41)
    
    Bumps [dawidd6/action-download-artifact](https://github.com/dawidd6/action-download-artifact) from 7 to 8.
    - [Release notes](https://github.com/dawidd6/action-download-artifact/releases)
    - [Commits](https://github.com/dawidd6/action-download-artifact/compare/v7...v8)
    
    ---
    updated-dependencies:
    - dependency-name: dawidd6/action-download-artifact
      dependency-type: direct:production
      update-type: version-update:semver-major
    ...
    
    Signed-off-by: dependabot[bot] <support@github.com>
    Co-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>

[33mcommit 4c1e24264edf26bccca27b250183287fe6af9578[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Sat Jan 25 11:16:37 2025 -0500

    add Bloch-Messiah decomposition (#37)

[33mcommit d17beaba9ae034f5860cdc3efc1128c3303932c4[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Sun Jan 19 19:02:41 2025 -0500

    space out docs (#36)

[33mcommit d0c314730bf5db948c954253193bbb989134b29a[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Tue Jan 14 21:44:13 2025 -0500

    Add symplectic analysis section to docs (#35)
    
    * update symplectic docs
    
    * update changelog and project.toml

[33mcommit c288e19024ee8c1aa3177f94b75f7d13930f9b3b[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Tue Jan 14 20:30:38 2025 -0500

    Add kwargs to `isapprox` on Gaussian types (#34)
    
    * add kwargs to `isapprox` on Gaussian types
    
    * add codecov
    
    * Add `changebasis` to public API (#33)
    
    * add changebasis to api
    
    * add more codecov
    
    * update project.toml and changelog
    
    * update changelog

[33mcommit bfafcd15e649495b7712083cac52e9685ce49131[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Tue Jan 14 20:16:28 2025 -0500

    Add `changebasis` to public API (#33)
    
    * add changebasis to api
    
    * add more codecov
    
    * update project.toml and changelog

[33mcommit f97533fc6ef9706355e007c3a048b82ff2d254dd[m
Author: github-actions[bot] <41898282+github-actions[bot]@users.noreply.github.com>
Date:   Tue Jan 14 16:32:32 2025 -0500

    CompatHelper: bump compat for Makie in [weakdeps] to 0.22, (keep existing compat) (#32)
    
    Co-authored-by: CompatHelper Julia <compathelper_noreply@julialang.org>

[33mcommit 9d438ba62ec6dff7d8af9b4a4654b54b0cf0b9fe[m
Author: apkille <killeandrew@gmail.com>
Date:   Thu Jan 9 18:19:27 2025 -0500

    update project.toml and changelog

[33mcommit 955d364264bf8dba0a8042e8e995f52bf6ab5045[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Thu Jan 9 18:16:43 2025 -0500

    Add symplectic polar decomposition (#28)
    
    * add polar decomp
    
    * add polar tests
    
    * export `polar` and fix tests

[33mcommit c5563f276745e458902bc1c612b387c8dd13dcc2[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Thu Jan 9 11:11:50 2025 -0500

    Add SymplecticFactorizations.jl as a dependency and williamson factorization support (#27)
    
    * williamson additions
    
    * update manifest

[33mcommit a1204298b3fab105319ada2d4ab5705b25fb56dc[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Mon Dec 30 17:40:13 2024 -0500

    Decrease max `nmode` size in tests to 5 to further prevent numerical instabilities (#26)
    
    * decrease max nmode size in tests to 5
    
    * update changelog

[33mcommit 03cbd2a200e9be85a8107956c9f24257fb423c18[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Mon Dec 30 17:13:45 2024 -0500

    fix `ptrace(::GaussianState, ::AbstractVector)` to preserve mode correlations (#25)
    
    * fix ptrace
    
    * update CHANGELOG and project.toml

[33mcommit d19112d426bb2232e4386154771edd3758865a33[m
Author: apkille <killeandrew@gmail.com>
Date:   Thu Dec 26 13:00:33 2024 -0500

    add example usage section to README

[33mcommit 938141b2092487d75af250f90d338663407a3273[m
Author: apkille <killeandrew@gmail.com>
Date:   Thu Dec 26 00:00:26 2024 -0500

    update project.toml and changelog to v1.2.3

[33mcommit 54875ae2f24d2c863a56a27717298594b709ade0[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Wed Dec 25 23:56:57 2024 -0500

    Propagate types in arrays to propagate ForwardDiff duals (#23)
    
    * promote types to propagate ForwardDiff duals
    
    * add ForwardDiff and FiniteDiff tests

[33mcommit ea87da6bfd96a4042ed6cc80a1558a7b687232c1[m
Author: apkille <killeandrew@gmail.com>
Date:   Wed Dec 25 04:27:30 2024 -0500

    reorganize manual in docs

[33mcommit c310e27e28965d7721ea8d769f317cf555dcb0ed[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Wed Dec 25 03:14:38 2024 -0500

    Add function for computing symplectic spectrum of a Gaussian state (#22)
    
    * symplectic spectrum
    
    * update project.toml and changelog.md
    
    * use det in test_states.jl
    
    * adjust tolerances for det test

[33mcommit 31864afe1084708c84d771723572861d6e867993[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Dec 21 23:39:32 2024 -0500

    trigger docs and rm sentence in docs

[33mcommit d96021bba52b46594c631ca751be46909aad04d2[m
Author: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>
Date:   Sat Dec 21 23:18:05 2024 -0500

    Bump dawidd6/action-download-artifact from 6 to 7 (#16)
    
    Bumps [dawidd6/action-download-artifact](https://github.com/dawidd6/action-download-artifact) from 6 to 7.
    - [Release notes](https://github.com/dawidd6/action-download-artifact/releases)
    - [Commits](https://github.com/dawidd6/action-download-artifact/compare/v6...v7)
    
    ---
    updated-dependencies:
    - dependency-name: dawidd6/action-download-artifact
      dependency-type: direct:production
      update-type: version-update:semver-major
    ...
    
    Signed-off-by: dependabot[bot] <support@github.com>
    Co-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>

[33mcommit 3d458081512b831a30b78ea78ee7a63b807f37ce[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Dec 17 20:04:54 2024 -0500

    update project.toml

[33mcommit 435e8fea6ae55ed2dd1e920120e5af252cd084f4[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Tue Dec 17 19:57:27 2024 -0500

    Implement `isgaussian` and `issymplectic` in tests and fix random generation of symplectic matrices for `QuadBlockBasis` (#21)
    
    * organize `test_random.jl` with new checks
    
    * change symplectic eig generation for block basis
    
    * update changelog

[33mcommit 53ee3429becd0766760a96dc6fd42e023b0abd46[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Tue Dec 17 18:54:29 2024 -0500

    Add `isgaussian` and `issymplectic` checks to public API (#20)
    
    * add `isgaussian` and `issymplectic`
    
    * add tests
    
    * reduce allocs
    
    * up project.toml and changelog
    
    * fix docstrings

[33mcommit 2a54aba9bd7462e39e2c21e256f9c622e208713f[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Dec 17 11:13:36 2024 -0500

    rm sentence in README

[33mcommit cc515c72c1d914087d36f8b9d38e0f629874fa43[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Dec 16 19:37:46 2024 -0500

    update project.toml

[33mcommit be28f241f3e610be24be0a81d1beb44abb06403a[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Mon Dec 16 19:35:34 2024 -0500

    Further establish symplectic interface (#19)
    
    * introduce `directsum` to API
    
    * add documentation section to manual
    
    * alter docs slightly
    
    * clean up symplectic basis table

[33mcommit 1d4d1cad1c93693d45ff900a46e352be8d97f7a1[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Dec 16 00:15:13 2024 -0500

    fix latex formatting error

[33mcommit f3f55cca8e7c599c72537537f2d14301fcaf0ea7[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Dec 14 08:47:27 2024 -0500

    print symplectic objects nicer

[33mcommit 36acb0f99f50a2724d8ed4791a6dc205a0c61a7d[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Sat Dec 14 07:53:29 2024 -0500

    Add symplectic representation feature (#15)
    
    * add `CanonicalForm` repr
    
    * blockform for symplecticform
    
    * use `basis` conventions from QuantumInterface
    
    * no longer make SymplecticBasis as subtype of Basis
    
    * update docstrings
    
    * add QuadBlockBasis features
    
    * fix tensor issues
    
    * more tests for QuadBlockBasis
    
    * update benchmark
    
    * add note on docs, will polish in future PR
    
    * rm loop errors in `twosqueeze` and `beamsplitter`
    
    * codecov
    
    * more codecov w/ states
    
    * more codecov!
    
    * this additional codecov will do...hopefully
    
    * update documentation
    
    * update CHANGELOG

[33mcommit b78da6eb83938d1fa11a5fe9594854640a68b5dc[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Mon Nov 18 23:11:46 2024 -0500

    Fix typo test error in Github workflow (#14)
    
    * fix typo test error in workflow
    
    * fix typos
    
    * fix benchmark.yml `judge` command
    
    * fix benchmark typo
    
    * change indices in `ptrace` benchmark
    
    * resolve manifest.toml
    
    * change dim to nmodes

[33mcommit 5033d731c3728d08d4abcf92d4a50266f9f7cbf6[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Mon Nov 18 22:07:08 2024 -0500

    Add benchmark suite to workflow (#13)
    
    * add benchmark suite
    
    * update CHANGELOG and project.toml
    
    * project.toml update

[33mcommit 88af062d08567dc852f55bd4c7bcc5e36e0ebfde[m
Author: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>
Date:   Mon Nov 18 14:49:40 2024 -0500

    Bump codecov/codecov-action from 4 to 5 (#12)
    
    Bumps [codecov/codecov-action](https://github.com/codecov/codecov-action) from 4 to 5.
    - [Release notes](https://github.com/codecov/codecov-action/releases)
    - [Changelog](https://github.com/codecov/codecov-action/blob/main/CHANGELOG.md)
    - [Commits](https://github.com/codecov/codecov-action/compare/v4...v5)
    
    ---
    updated-dependencies:
    - dependency-name: codecov/codecov-action
      dependency-type: direct:production
      update-type: version-update:semver-major
    ...
    
    Signed-off-by: dependabot[bot] <support@github.com>
    Co-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>

[33mcommit d9580fcb2a0413166d7d6daa820e79d54b691da5[m
Author: apkille <killeandrew@gmail.com>
Date:   Sun Nov 17 21:28:14 2024 -0500

    add note about symplectic forms in docs

[33mcommit 6a02a85d0d0c202e773d831dfe4adfc95afce656[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Sun Nov 17 20:56:50 2024 -0500

    Generate valid Gaussian states, unitaries, and channels (#10)
    
    * add random features
    
    * add more tests
    
    * fix
    
    * change basis of symplecticform
    
    * update project.toml and CHANGELOG
    
    * alter Float32 test
    
    * add codecov for `randsymplectic`
    
    * change floating point tests for rand objects

[33mcommit 4d662d31227a9e2e636b802e345e0068ac9ff506[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Tue Nov 5 23:53:39 2024 -0500

    add positive-definite condition (#9)

[33mcommit f9ead0ef5d27fd03ac3a2ff3a371f9563bb22cab[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Nov 5 16:18:15 2024 -0500

    fix `GaussianChannel` docstring

[33mcommit feaff1d6243f0be339d42a83e8581164c698640f[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Nov 5 15:05:14 2024 -0500

    add docstrings to `prob` and `output`

[33mcommit fe61484408acb5b4e52d97c9506ef2dca978c03f[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Tue Nov 5 00:42:49 2024 -0500

    Add `prob` function for `Generaldyne` type (#8)
    
    * add `prob` function for `Generaldyne` type
    
    * update CHANGELOG.md and Project.toml

[33mcommit 762b99a13b85623d67203783d0a50779f5575fd5[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Nov 4 19:57:16 2024 -0500

    documentation small edits

[33mcommit 177a2510a0b842dbc37cadfbaee28a8776dbabb1[m
Author: apkille <killeandrew@gmail.com>
Date:   Sun Nov 3 20:53:49 2024 -0500

    update Project.toml

[33mcommit d52b94b7529dc9f1441e666c40be55b3441bfe96[m
Author: apkille <killeandrew@gmail.com>
Date:   Sun Nov 3 19:54:36 2024 -0500

    update CHANGELOG.md date

[33mcommit 378e147229eeb37f208ef81faf11782764db7a51[m
Author: apkille <killeandrew@gmail.com>
Date:   Sun Nov 3 19:47:35 2024 -0500

    update Project.toml

[33mcommit ea44c7554105810fe82aeee80f5dd7a544d68de5[m
Author: apkille <killeandrew@gmail.com>
Date:   Sun Nov 3 19:39:43 2024 -0500

    update CHANGELOG.md

[33mcommit ac440724d997223c7fd4eff935512a1a8bd8fc06[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Sun Nov 3 19:29:18 2024 -0500

    add'randstate' and 'randchannel' functions (#7)

[33mcommit 56679c0e88c80b0cef181463ac9ecc5ab88c50c2[m
Author: apkille <killeandrew@gmail.com>
Date:   Thu Oct 31 23:30:50 2024 -0400

    def `isapprox` for operators

[33mcommit ab75fd8661961a9ff6f6233c7410a142ccc06baa[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Thu Oct 31 23:22:18 2024 -0400

    Rename `outcome` function with `output` (#6)
    
    * replace 'outcome' with 'output'
    
    * update CHANGELOG

[33mcommit ade660354359d7761308d0bc9bc144c7792b8a29[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Tue Oct 29 22:24:05 2024 -0400

    Add Generaldyne measurements (#5)
    
    * initial steps
    
    * BlockArray additions
    
    * tests
    
    * add tests
    
    * update changelog

[33mcommit 42cdb2ad1ec20c7d5405a2bd661de787de04f3d4[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Oct 28 11:08:41 2024 -0400

    bump Julia compat in Project.toml

[33mcommit 9e0a0c66b6b1951aac0733c40c68353ab82dd273[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Oct 22 23:59:45 2024 -0400

    readings in readme

[33mcommit 5949248a44aa88556a269c22007073e75f52a725[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Oct 22 23:54:01 2024 -0400

    add install instructions

[33mcommit 9bcf5c137409f7c12d82548c3ac35169f7660521[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Oct 21 05:16:53 2024 -0400

    update note in tutorial

[33mcommit 0279cd1c4f9c5a946181b9c0b8abe16ef81149a2[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Oct 21 05:16:32 2024 -0400

    update note in tutorial

[33mcommit 722cbee14f22146fb0aa6a271174155ef84ba77a[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Oct 21 05:04:32 2024 -0400

    operator instances in readme

[33mcommit 0d539a859df5f8fc473f4c1e2e767638ec8235b3[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Oct 21 05:03:44 2024 -0400

    click me edit

[33mcommit aace7c821ce29f435e903ee0d2a7c7d801c17dc9[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Oct 21 05:02:01 2024 -0400

    add ex for Gaussian objs

[33mcommit 08f7255f54cd5d1b25ba45222818dc3667492e69[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Oct 21 04:16:16 2024 -0400

    update project.toml

[33mcommit 60d453d0e8b705a316d8da808de155aeaad64590[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Oct 21 04:05:08 2024 -0400

    add note about contributing and getting started

[33mcommit 5a12adccf0649f9a272ac24060492fa8e58cd774[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Oct 21 03:17:18 2024 -0400

    add ex in custom array tutorial

[33mcommit d1f5770072b1eda13c6dc0c3a444d5b51a21e08a[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Oct 19 15:00:13 2024 -0400

    add docs for Custom Array types

[33mcommit e4bca87f2b99f6cac3756089c381eeeadba18dd5[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Oct 19 11:38:18 2024 -0400

    polish Makie tutorial with another example

[33mcommit 5b031043bff833ee78f5d2729c9e2823b4e4182a[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Sat Oct 19 10:47:33 2024 -0400

    Remove StaticArrays as dependency and add as extension (#4)
    
    * create StaticArrays ext
    
    * update changelog and project.toml
    
    * rm uncovered methods in StaticArraysExt/utils.jl

[33mcommit 28227261c36790cb3c84149a7d2a3c20adf2e187[m
Author: github-actions[bot] <41898282+github-actions[bot]@users.noreply.github.com>
Date:   Thu Oct 17 22:32:38 2024 -0400

    CompatHelper: add new compat entry for Makie in [weakdeps] at version 0.21, (keep existing compat) (#3)
    
    Co-authored-by: CompatHelper Julia <compathelper_noreply@julialang.org>

[33mcommit 63ad3b9b7a4225188d1dd2dd0695590fb331120f[m
Author: apkille <killeandrew@gmail.com>
Date:   Thu Oct 17 19:05:08 2024 -0400

    add more test cov for channels

[33mcommit 22cac16c9b53661b130608dd45bcb4bd6d768d58[m
Author: apkille <killeandrew@gmail.com>
Date:   Wed Oct 16 22:15:40 2024 -0400

    bump version number of Project.toml

[33mcommit 1c725fefc2b85be9f23c0ad7cd1bdeb48e9b9f71[m
Author: apkille <killeandrew@gmail.com>
Date:   Wed Oct 16 22:01:10 2024 -0400

    rm Makie.jl from depend and polish Makie docs

[33mcommit 50baf235dc86b1fe0eae1071a633022eb992a9c6[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Oct 12 17:40:53 2024 -0400

    add CairoMakie to docs Project.toml

[33mcommit 7860ed5f540e1b71dff2a47223ed5e7b03b186ad[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Thu Oct 10 00:25:11 2024 -0400

    add Makie heatmap attributes (#2)
    
    * add Makie heatmap attributes
    
    * add tutorial for Makie
    
    * add phase space coord comment
    
    * update changelog

[33mcommit a85d1e0e7e2dd71a64e8c92d26172636cc36009b[m
Author: apkille <killeandrew@gmail.com>
Date:   Thu Sep 26 19:28:22 2024 -0400

    reduce allocs in two state predef objects

[33mcommit 0d0d1bb9f6a267db46682253230e53eb41c4e5f3[m
Author: apkille <killeandrew@gmail.com>
Date:   Thu Sep 26 15:56:59 2024 -0400

    add Gabs.jl  pkg name in installation

[33mcommit c1757473019c7274d077af57f0e018e90e4ec262[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Sep 24 14:01:06 2024 -0400

    add CHANGELOG for first release

[33mcommit a8c8321a57b4028fdd92c473578b286ff7c7c738[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Sep 21 23:23:18 2024 -0400

    add learning resources in docs

[33mcommit 81c6cfd35388452104f13e7543abf4c1ca6d2aff[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Sep 21 13:52:59 2024 -0400

    improve README clarity

[33mcommit 326ae206908d48a0b6d473d90fe2ecaf1a118ac4[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Sep 21 12:35:53 2024 -0400

    include README description

[33mcommit 72643f1a0edf44eda08235f107d714920514ef79[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Sep 21 11:55:16 2024 -0400

    trigger docs

[33mcommit 57e8b8e42048837776e56332cf043bc67403b5c7[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Sep 21 00:40:50 2024 -0400

    trigger doc build

[33mcommit 85d492cfd71d871ef40f0ea8d764924a1c29ad44[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Sep 21 00:28:55 2024 -0400

    add .git to deploydocs repo

[33mcommit 922a87a88f014003943516b209d46ed313c05811[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Sep 21 00:11:15 2024 -0400

    set deploydocs to main branch

[33mcommit 7c38f3d42d7e2515120624288a3bc3c558e93c66[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Sep 21 00:02:47 2024 -0400

    add CairoMakie to test project.toml

[33mcommit 3c1bd9dc5d2e2d41803667afc4302498b6f2b4d4[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 20 23:59:18 2024 -0400

    add Documenter to test Project.toml

[33mcommit c98aa44dd11b4f0a1a071448b9d9d7dadf227189[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 20 23:56:40 2024 -0400

    rm manifest.toml

[33mcommit 25d9cf9c3283b0347a4320e27fc217c1d419691c[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 20 23:42:42 2024 -0400

    canonical format for docs

[33mcommit 6c244e39a4406807d21e455e15cb8801a1445e1d[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 20 23:39:56 2024 -0400

    deploy docs

[33mcommit 478892ee27728c9fcd91d01d91ab1bd92fc1b72a[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 20 22:49:46 2024 -0400

    update to version 1.0.0

[33mcommit fb78db199714702220a5a33fa94f4e1d33b05568[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 20 22:34:29 2024 -0400

    add badges

[33mcommit 372de5f714a614e7474f327583daba658012a8eb[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 20 20:53:36 2024 -0400

    rm doctest fix

[33mcommit b30314fc82f8cc3a9e37c53423f22dc94bc4cf6f[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 20 20:52:53 2024 -0400

    improve Base.summary for Gaussian types

[33mcommit c660c9de84990ab5cac41ab90b9a820267045d3e[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 20 20:15:13 2024 -0400

    improve `ptrace` and `tensor` w/ nmodes

[33mcommit cfce324d2236225409d4fd114ab2ffe388fc3ae2[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 20 19:48:10 2024 -0400

    add nmodes field to Gaussian types to reduce alloc

[33mcommit 3bfc8d6ed6bbf427226c09e1b270cb848644edc9[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 20 18:41:28 2024 -0400

    make ptrace preserve array types

[33mcommit 222e5f73efd29d1a0b1210515432aea2df6a9dd6[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Sep 17 22:27:17 2024 -0400

    open up `tensor` to custom array types

[33mcommit f3744b5732c114f23b10117e2c37e0da3a19f6ea[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Sep 17 16:24:04 2024 -0400

    rm apply in place for Base.:(*)

[33mcommit 8756fee81616737be4920324d534ae69a14a3df7[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Sep 17 09:15:20 2024 -0400

    save states.jl

[33mcommit 03350636e024ee45e390c95cb645000579450bfa[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Sep 17 09:14:26 2024 -0400

    direct sum change to tensor for clarity

[33mcommit 320d448667540090ed66a55e75a6f2efcb8bc7fe[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 13 21:06:07 2024 -0400

    makie convert typo fix

[33mcommit 5c61c8d6d46c800284efa2970eab94041ce19e88[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 13 06:05:49 2024 -0400

    add names for tutorials section

[33mcommit 628aec2022950b80a7e7411b9c08cb88ea5bdc47[m
Author: apkille <killeandrew@gmail.com>
Date:   Thu Sep 12 13:52:53 2024 -0400

    add predef channels to manual

[33mcommit 4860b216b656f2e3c552e1552ef22ce21ffb6b2f[m
Author: apkille <killeandrew@gmail.com>
Date:   Wed Sep 11 21:21:26 2024 -0400

    add amplifier channel

[33mcommit 43deebb9311ba6b10ad9ed2031651964970ddd6e[m
Author: apkille <killeandrew@gmail.com>
Date:   Wed Sep 11 21:11:54 2024 -0400

    attenuator docstring

[33mcommit 446e3ee75153ce9e5daa9b5124683b944f24d5b1[m
Author: apkille <killeandrew@gmail.com>
Date:   Wed Sep 11 21:11:37 2024 -0400

    better attenuator docstring

[33mcommit 10770350508936083767162765d428094cd36f6e[m
Author: apkille <killeandrew@gmail.com>
Date:   Wed Sep 11 21:04:02 2024 -0400

    add attenuator channel

[33mcommit 9240774abb765149786b03f5f2d2554830a387dd[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Sep 10 20:28:05 2024 -0400

    fixing directsum docstring

[33mcommit 0a554c9bf67b9b066c6d4d83a0c608587bd7300e[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Sep 10 19:17:44 2024 -0400

    implement ptrace

[33mcommit 380f4993356ecc7ac41d79dab83a537d08cc6b4c[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Sep 10 11:26:03 2024 -0400

    define direct sum util function

[33mcommit 35b4abc6907fdf0f3d05fbd24db00d5e1974ee13[m
Author: apkille <killeandrew@gmail.com>
Date:   Sat Sep 7 08:24:19 2024 -0400

    using reorganization in tests

[33mcommit 79308e24afd76b2ac11634e336b4156ce649b3d4[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 6 21:48:23 2024 -0400

    fix doi link in references

[33mcommit e6912f2d3c0902d6f9ed6ecb7f3191351b3e85fc[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 6 21:21:59 2024 -0400

    update GaussianChannel part of manual

[33mcommit be0d52a763ca392f99570a7d5898aeac9af42d72[m
Author: apkille <killeandrew@gmail.com>
Date:   Fri Sep 6 19:37:45 2024 -0400

    rearrange operations on Gaussian objects

[33mcommit f7556ef875c7dc0677740b4734089a0a63618f2a[m
Author: apkille <killeandrew@gmail.com>
Date:   Thu Sep 5 17:55:36 2024 -0400

    add makie error and test

[33mcommit 6d7cd8a9dc1d66f63632b164ee9dd30cb10cab81[m
Author: apkille <killeandrew@gmail.com>
Date:   Thu Sep 5 13:31:39 2024 -0400

    add Makie ext

[33mcommit f511cd2a2f8a52a1b44b6cae3a3c397cc2fc6b02[m
Author: apkille <killeandrew@gmail.com>
Date:   Wed Sep 4 21:56:40 2024 -0400

    include citations for review

[33mcommit 36d11ab0e2093ba3cdc7a37133720eb6af9fd31e[m
Author: apkille <killeandrew@gmail.com>
Date:   Wed Sep 4 21:53:52 2024 -0400

    wigner functions

[33mcommit bcece3eda3cad8e969386a439a11063f74ab0639[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Sep 3 14:31:29 2024 -0400

    add getting started and tutorial pages

[33mcommit ace586ecae5630b7679a9e82fc9c4337ff4a15bb[m
Author: apkille <killeandrew@gmail.com>
Date:   Tue Sep 3 12:39:41 2024 -0400

    docs for Gaussian Unitarites

[33mcommit 90f530ad65262fcbaeca3282f37bcb84ad5e71ae[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Sep 2 20:18:25 2024 -0400

    create low-level manual and zoo

[33mcommit 049654e90d40c933c3ddbe4b049aba14d138d4d0[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Sep 2 17:26:05 2024 -0400

    add docstrings with tests

[33mcommit c127a75b8588ebfcba22a69c0f71bcfad36e6006[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Sep 2 14:28:12 2024 -0400

    add GaussianUnitary

[33mcommit 7f34a5f69d28ed8a92d1f01a15cedce19295db34[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Sep 2 11:30:16 2024 -0400

    add docs folder

[33mcommit fa5a46e4a169bb448b78d7c7c1ae44699f95b2b7[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Sep 2 11:08:07 2024 -0400

    add github workflows

[33mcommit 7dd2124454b5a899806ea1c00ef6b5168cb7bc73[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Sep 2 11:02:05 2024 -0400

    add ci buildkite

[33mcommit bd952e06a7e1cb013d431f8452e4528381b8cb71[m
Author: apkille <killeandrew@gmail.com>
Date:   Mon Sep 2 04:45:50 2024 -0400

    add Gaussian states and operators

[33mcommit 06c2821a46c930e985848cbb162c00a1220e6519[m
Author: apkille <killeandrew@gmail.com>
Date:   Sun Sep 1 20:27:53 2024 -0400

    create pkg skeleton

[33mcommit 65cdf5ced8ab389713011b2be1d780cbaac80846[m
Author: Andrew Kille <68079167+apkille@users.noreply.github.com>
Date:   Sun Sep 1 20:17:36 2024 -0400

    Initial commit
