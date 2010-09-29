#called when package is loaded :
.onLoad = function (lib, pkg) {
    library.dynam(pkgnm(), pkg, lib)
}

#called when namespace is unloaded ( unloadNamespace("package:pkg_name") )
.onUnload = function(path) {
    library.dynam.unload(pkgnm(), path)
}

#called when package is detached ( detach("package:pkg_name") )
.Last.lib = function(path) {
    library.dynam.unload(pkgnm(), path)
}

