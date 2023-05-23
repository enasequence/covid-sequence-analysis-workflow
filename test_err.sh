err_file=log.err
echo "hello" 2>> ${err_file}
echo "hello2" 2>> ${err_file}

if [ ! -s ${err_file} ]; then
  rm ${err_file}
fi
command_name 2>> ${err_file}
nosuchcommand1 2>> ${err_file}
nosuchcommand2 2>> ${err_file}