# This file has utilities to download DFMO and TESTMO from github.
# We make some assumptions on the structure of both packages so that the functions
# are relatively simple and we do not have to care about edge cases etc.

function download_testmo(testmo_path; overwrite=false)
    return download_lib(testmo_path, TESTMO_TARBALL_URL; overwrite)
end

function download_dfmo(dfmo_path; overwrite=false)
    return download_lib(dfmo_path, DFMO_TARBALL_URL; overwrite)
end

function download_lib(lib_path, lib_url; overwrite=false)
    ## make sure there is an empty target directory:
    _lib_path = prepare_target_path(lib_path; overwrite)
    isnothing(_lib_path) && return _lib_path
    
    ## download and extract compressed tarball
    _lib_path = download_url_to_path(lib_url, _lib_path)
    isnothing(_lib_path) && return _lib_path

    ## now the directory has structure `_lib_path/DFMO-main` or similar
    ## we move all files of the single subfolder to its parent `_lib_path`:
    contents = readdir(_lib_path)
    dl_dir = joinpath(_lib_path, only(contents))
    for fn in readdir(dl_dir)
        mv(joinpath(dl_dir, fn), joinpath(_lib_path, fn); force=true)
    end
    ## remove empty subfolder
    rm(dl_dir)
    ## return 
    return _lib_path
end

function download_url_to_path(url, full_trgt_path = tempname())
    if isdir(full_trgt_path) && isempty(readdir(full_trgt_path))
        io = IOBuffer()
        DL.download(url, io)
        seekstart(io)
        tar = GzipDecompressorStream(io)
        Tar.extract(tar, full_trgt_path)
        close(io)
        close(tar)
        
        return full_trgt_path
    else
        @warn "Download canceled -- Target `$full_trgt_path` either not existent or not empty."
        return nothing
    end
end