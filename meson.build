project('otfp', 'c',
     version: run_command(meson.current_source_dir()/'version.sh', check: true).stdout().strip(),
     meson_version : '>=0.57.0',
)

meson.add_dist_script('version.sh')

subdir('src')

bash = find_program('bash', required: false)
if bash.found()
  test(
    'examples',
    bash,
    args: [meson.current_source_dir() / 'examples' / 'test.sh'],
    workdir: meson.current_source_dir() / 'examples',
    timeout: 300,
    is_parallel: false,
  )
endif
