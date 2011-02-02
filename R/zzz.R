#called when package is detached ( detach("package:pkg_name") )
.Last.lib = function(path) {
    library.dynam.unload(pkgnm(), path)
}

