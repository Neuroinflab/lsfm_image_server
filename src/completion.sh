# bash completion support for ./lmdf_cli.py

_complete-lmdf_clipy()
{
  local start cur prev opts
  COMPREPLY=()
  start="${COMP_WORDS[@]:0:1}"
  if [ ${COMP_CWORD} -gt 1  ]  ; then
    start="${COMP_WORDS[@]:0:2}"
  fi
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"

  opts=""

  if [[ "$prev" == "--hdf-path" ]] ; then
    COMPREPLY=( $(compgen -o plusdirs -o nospace -f -X '!*.h5' -- ${cur}) )
    return 0
  fi

  if [[ "$prev" == "--affine-path" ]] ; then
    COMPREPLY=( $(compgen -o plusdirs -o nospace -f -X '!*.txt' -- ${cur}) )
    return 0
  fi

  if [[ "$prev" == "--metadata-path" ]] ; then
    COMPREPLY=( $(compgen -o plusdirs -o nospace -f -X '!*.json' -- ${cur}) )
    return 0
  fi

  if [[ "$prev" == "--image-path" ]] ; then
    COMPREPLY=( $(compgen -o plusdirs -o nospace -f -X '!*.tif' -- ${cur}) )
    return 0
  fi

  if [[ "$start" == "./lmdf_cli.py info" ]] ; then
    opts="--channel-name --hdf-path"
  fi

  if [[ "$start" == "./lmdf-cli.py add-metadata" ]] ; then
    opts="--hdf-path --m-key --m-value"
  fi

  if [[ "$start" == "./lmdf-cli.py export" ]] ; then
    opts="--channel-name --grid-size --hdf-path --input-orientation --input-resolution-level --output-path --output-resolution --overlap-mm --phys-origin --phys-size --region-id --segmentation-name"
  fi

  if [[ "$start" == "./lmdf_cli.py export" ]] ; then
    opts="--channel-name --grid-size --hdf-path --input-orientation --input-resolution-level --output-path --output-resolution --overlap-mm --phys-origin --phys-size --region-id --segmentation-name"
  fi

  if [[ "$start" == "./lmdf-cli.py add-transform" ]] ; then
    opts="--invert --name --ordn --transform-type"
  fi

  if [[ "$start" == "./lmdf-cli.py info" ]] ; then
    opts="--channel-name --hdf-path"
  fi

  if [[ "$start" == "./lmdf_cli.py add-metadata" ]] ; then
    opts="--hdf-path --m-key --m-value"
  fi

  if [[ "$start" == "./lmdf-cli.py export-slices" ]] ; then
    opts="--axis --channel-name --extract-roi --hdf-path --input-orientation --input-resolution-level --output-path --slicing-range"
  fi

  if [[ "$start" == "./lmdf-cli.py write-affine" ]] ; then
    opts="--affine-name --affine-path --channel-name --hdf-path"
  fi

  if [[ "$start" == "./lmdf-cli.py write" ]] ; then
    opts="--bdv-xml --channel-name --file-name-format --hdf-path --image-path --is-multichannel --is-segmentation --metadata-path --slab-memory-size"
  fi

  if [[ "$start" == "./lmdf_cli.py show-log" ]] ; then
    opts="--hdf-path"
  fi

  if [[ "$start" == "./lmdf_cli.py export-slices" ]] ; then
    opts="--axis --channel-name --extract-roi --hdf-path --input-orientation --input-resolution-level --output-path --slicing-range"
  fi

  if [[ "$start" == "./lmdf-cli.py show-log" ]] ; then
    opts="--hdf-path"
  fi

  if [[ "$start" == "./lmdf_cli.py" ]] ; then
    opts="add-metadata add-transform export export-slices info show-log write write-affine"
  fi

  if [[ "$start" == "./lmdf-cli.py" ]] ; then
    opts="add-metadata add-transform export export-slices info show-log write write-affine"
  fi

  if [[ "$start" == "./lmdf_cli.py write" ]] ; then
    opts="--bdv-xml --channel-name --file-name-format --hdf-path --image-path --is-multichannel --is-segmentation --metadata-path --slab-memory-size"
  fi

  if [[ "$start" == "./lmdf_cli.py add-transform" ]] ; then
    opts="--invert --name --ordn --transform-type"
  fi

  if [[ "$start" == "./lmdf_cli.py write-affine" ]] ; then
    opts="--affine-name --affine-path --channel-name --hdf-path"
  fi

  COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
  return 0
}

complete -F _complete-lmdf_clipy ./lmdf_cli.py

