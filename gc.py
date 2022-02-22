import os

from google.cloud import storage


def get_paths(client, bucket_dir):
    all_paths = dict()

    for blob in client.list_blobs(bucket_dir):
        folder = os.path.dirname(blob.name)
        if 'illumina' not in folder:
            continue
        else:
            folder_spl = folder.split('/')
            key_dir = folder_spl[0]

            if len(folder_spl) in [1, 2, 3] and key_dir not in all_paths:
                all_paths[key_dir] = {'tar_count': 0, 'folder_files_count': 0, 'folder_files': dict()}
                print(f"{key_dir} created in dict by folder {folder}")

            if len(folder_spl) == 1:
                continue
            elif len(folder_spl) == 2:
                if 'results' not in folder:
                    print(f"'results' not in {folder}")
                else:
                    if 'tar' in blob.name:
                        all_paths[key_dir]['tar_count'] += 1
                    else:
                        print(f"current folder is {blob.name}")
            elif len(folder_spl) == 3:
                if 'results' not in folder:
                    print(f"'results' not in {folder}")
                else:
                    all_paths[key_dir]['folder_files_count'] += 1
                    if folder_spl[2] not in all_paths[key_dir]['folder_files']:
                        all_paths[key_dir]['folder_files'][folder_spl[2]] = 0
                    all_paths[key_dir]['folder_files'][folder_spl[2]] += 1
            else:
                print(f"err with path length: {folder}")

    return all_paths


def get_prev_paths(client, bucket_dir):
    all_prev_paths = dict()

    for blob in client.list_blobs(bucket_dir, prefix='2021-12-13'):
        folder = os.path.dirname(blob.name)
        if 'illumina' not in folder:
            continue
        else:
            folder_spl = folder.split('/')
            key_dir = folder_spl[1]

            if len(folder_spl) == 3:
                if 'results' not in folder:
                    print(f"'results' not in {folder}")
                else:
                    if 'tar' in blob.name:
                        if key_dir not in all_prev_paths:
                            all_prev_paths[key_dir] = {'tar_count': 0, 'folder_files_count': 0, 'folder_files': dict()}
                        all_prev_paths[key_dir]['tar_count'] += 1
                    else:
                        print(f"current folder is {blob.name}")
            elif len(folder_spl) == 4:
                if 'results' not in folder:
                    print(f"'results' not in {folder}")
                else:
                    all_prev_paths[key_dir]['folder_files_count'] += 1
                    if folder_spl[3] not in all_prev_paths[key_dir]['folder_files']:
                        all_prev_paths[key_dir]['folder_files'][folder_spl[3]] = 0
                    all_prev_paths[key_dir]['folder_files'][folder_spl[3]] += 1
            else:
                print(f"err with path length: {folder}")

    return all_prev_paths


if __name__ == "__main__":
    client = storage.Client()

    bucket_dir = '2022-01-25-upd-bucket'
    all_paths = get_paths(client=client, bucket_dir=bucket_dir)

    print(all_paths.keys())
    print(len(all_paths))
    for key in all_paths:
        print(key, all_paths[key]['tar_count'], all_paths[key]['folder_files_count'] / 2)
        '''if all_paths[key]['tar_count'] != all_paths[key]['folder_files_count'] / 2:
            dict_test = all_paths[key]['folder_files']
            for k, v in dict_test.items():
                if v != 2:
                    print(k, v)'''
    print('---------')

    bucket_prev_dir = 'prj-int-dev-covid19-nf-gls'
    all_prev_paths = get_prev_paths(client=client, bucket_dir=bucket_prev_dir)

    print(all_prev_paths.keys())
    print(len(all_prev_paths))
    for key in all_prev_paths:
        print(key, all_prev_paths[key]['tar_count'], all_prev_paths[key]['folder_files_count'] / 6)
        '''if all_prev_paths[key]['tar_count'] != all_prev_paths[key]['folder_files_count'] / 6:
            dict_test = all_prev_paths[key]['folder_files']
            for k, v in dict_test.items():
                if v != 6:
                    print(k, v)'''
    print('---------')

    for key in ['illumina_19', 'illumina_20', 'illumina_21']:
        keys_curr = set(all_paths[key]['folder_files'].keys())
        keys_prev = set(all_prev_paths[key]['folder_files'].keys())
        print(key, keys_curr-keys_prev, keys_prev-keys_curr)
