#!/usr/bin/env python3

# 
# This installs the module dune-wrapper and its dependencies.
# The exact revisions used are listed in the table below.
# However, note that this script may also apply further patches.
# If so, all patches are required to be the current folder, or,
# in the one that you specified as argument to this script.
# 
# 
# |      module name      |      branch name      |                 commit sha                 |         commit date         |
# |-----------------------|-----------------------|--------------------------------------------|-----------------------------|
# |      dune-spgrid      |  origin/releases/2.9  |  8f7e156e8bb68543e195e44a012bc8453e18ddac  |  2022-11-05 10:40:40 +0000  |
# |  dune-localfunctions  |  origin/releases/2.9  |  fc5c8050452c59c335e6be28afb761bdbd4479ad  |  2023-08-21 22:09:53 +0000  |
# |         dumux         |  origin/releases/3.8  |  c8f61c1f81ca511415c656e834cc0ded17572025  |  2023-12-01 10:12:26 +0000  |
# |      dune-common      |  origin/releases/2.9  |  51f8c97a4fd06bbd60164a6d241db2a605d3e139  |  2023-10-31 14:56:11 +0000  |
# |     dumux-precice     |     origin/develop    |  336b7bb2c1419332cd138a1a05666a82166d0fc6  |  2023-12-11 04:10:11 +0000  |
# |       dune-grid       |  origin/releases/2.9  |  cd3ac9e22402301a8c8290aa721e30af6d4b312d  |  2023-06-21 11:24:38 +0000  |
# |       dune-istl       |  origin/releases/2.9  |  9cb74375ff0970384c99a7db83ce5f292d06b1a4  |  2023-08-21 22:06:08 +0000  |
# |      dune-subgrid     |     origin/master     |  41ab447c59ea508c4b965be935b81928e7985a6b  |  2022-09-25 23:18:45 +0000  |
# |    dumux-phasefield   |  origin/cell_problems |  47ad96f8c4924853bd8228a434bf1713b4a848b1  |  2023-04-14 08:09:45 +0000  |
# |     dune-geometry     |  origin/releases/2.9  |  1051f4d4e7ea10ec7787e49974c4ff5d4debc176  |  2023-07-15 07:21:03 +0000  |

import os
import sys
import subprocess

top = "."
os.makedirs(top, exist_ok=True)


def runFromSubFolder(cmd, subFolder):
    folder = os.path.join(top, subFolder)
    try:
        subprocess.run(cmd, cwd=folder, check=True)
    except Exception as e:
        cmdString = ' '.join(cmd)
        sys.exit(
            "Error when calling:\n{}\n-> folder: {}\n-> error: {}"
            .format(cmdString, folder, str(e))
        )


def installModule(subFolder, url, branch, revision):
    targetFolder = url.split("/")[-1]
    if targetFolder.endswith(".git"):
        targetFolder = targetFolder[:-4]
    if not os.path.exists(targetFolder):
        runFromSubFolder(['git', 'clone', url, targetFolder], '.')
        runFromSubFolder(['git', 'checkout', branch], subFolder)
        runFromSubFolder(['git', 'reset', '--hard', revision], subFolder)
    else:
        print(
            f"Skip cloning {url} since target '{targetFolder}' already exists."
        )


def applyPatch(subFolder, patch):
    sfPath = os.path.join(top, subFolder)
    patchPath = os.path.join(sfPath, 'tmp.patch')
    with open(patchPath, 'w') as patchFile:
        patchFile.write(patch)
    runFromSubFolder(['git', 'apply', 'tmp.patch'], subFolder)
    os.remove(patchPath)

print("Installing dune-spgrid")
installModule("dune-spgrid", "https://gitlab.dune-project.org/extensions/dune-spgrid.git", "origin/releases/2.9", "8f7e156e8bb68543e195e44a012bc8453e18ddac", )

print("Installing dune-localfunctions")
installModule("dune-localfunctions", "https://gitlab.dune-project.org/core/dune-localfunctions.git", "origin/releases/2.9", "fc5c8050452c59c335e6be28afb761bdbd4479ad", )

print("Installing dumux")
installModule("dumux", "git@git.iws.uni-stuttgart.de:dumux-repositories/dumux.git", "origin/releases/3.8", "c8f61c1f81ca511415c656e834cc0ded17572025", )

print("Installing dune-common")
installModule("dune-common", "https://gitlab.dune-project.org/core/dune-common.git", "origin/releases/2.9", "51f8c97a4fd06bbd60164a6d241db2a605d3e139", )

print("Installing dumux-precice")
installModule("dumux-adapter", "https://github.com/precice/dumux-adapter.git", "origin/develop", "336b7bb2c1419332cd138a1a05666a82166d0fc6", )

print("Installing dune-grid")
installModule("dune-grid", "https://gitlab.dune-project.org/core/dune-grid.git", "origin/releases/2.9", "cd3ac9e22402301a8c8290aa721e30af6d4b312d", )

print("Installing dune-istl")
installModule("dune-istl", "https://gitlab.dune-project.org/core/dune-istl.git", "origin/releases/2.9", "9cb74375ff0970384c99a7db83ce5f292d06b1a4", )

print("Installing dune-subgrid")
installModule("dune-subgrid", "ssh://git@gitlab.dune-project.org:22022/extensions/dune-subgrid.git", "origin/master", "41ab447c59ea508c4b965be935b81928e7985a6b", )

print("Installing dune-geometry")
installModule("dune-geometry", "https://gitlab.dune-project.org/core/dune-geometry.git", "origin/releases/2.9", "1051f4d4e7ea10ec7787e49974c4ff5d4debc176", )

print("Installing dumux-phasefield")
installModule("dumux-phasefield", "git@git.iws.uni-stuttgart.de:dumux-appl/dumux-phasefield.git", "origin/cell_problems", "47ad96f8c4924853bd8228a434bf1713b4a848b1", )

print("Applying patch for unpublished commits in dune-subgrid")
patch = """
From 6d01ec3c1e8656bfdf0c1339c86445a4be0c8ca0 Mon Sep 17 00:00:00 2001
From: Ned Coltman <edward.coltman@iws.uni-stuttgart.de>
Date: Wed, 25 May 2022 16:03:59 +0200
Subject: [PATCH] [TEMP] This works for what I need it to do, but the YASP and
 ALU tests fail

---
 dune/subgrid/subgrid/subgridintersection.hh | 20 ++++++++++++++++----
 1 file changed, 16 insertions(+), 4 deletions(-)

diff --git a/dune/subgrid/subgrid/subgridintersection.hh b/dune/subgrid/subgrid/subgridintersection.hh
index 8e231df..d48b643 100644
--- a/dune/subgrid/subgrid/subgridintersection.hh
+++ b/dune/subgrid/subgrid/subgridintersection.hh
@@ -167,18 +167,30 @@ public:
     //! Return outside element, throws exception for boundary elements
     Entity outside() const {
         if (!neighbor())
-            DUNE_THROW(GridError, \"There is no neighbor!\");
-        return outside_;
+            DUNE_THROW(GridError,\"no in neighbor found in outside()\");
+
+        HostEntity outside = insideIntersect_.outside();
+        if (!this->subGrid_->template contains<0>(outside))
+            DUNE_THROW(GridError,\"no neighbor found in outside()\");
+
+        return SubEntity(this->subGrid_, std::move(outside));
     }

     //! return true if across the edge a neighbor exists
     bool neighbor () const {
-        return inside_ != outside_;
+        if (insideIntersect_.neighbor())
+            return (this->subGrid_->template contains<0>(insideIntersect_.outside()));
+        else
+            return false;
     }

     //! return true if intersection is with boundary.
     bool boundary () const {
-        return !neighbor();
+        if (neighbor())
+            return ( insideIntersect_.boundary()
+                  || !this->subGrid_->template contains<0>(insideIntersect_.outside()));
+        else
+            return true;
     }

     //! local number of codim 1 entity in neighbor where intersection is contained
--
2.34.1
"""
applyPatch("dune-subgrid", patch)

print("Configuring project")
runFromSubFolder(
    ['./dune-common/bin/dunecontrol', '--opts=dune-wrapper/cmake.opts', 'all'],
    '.'
)
