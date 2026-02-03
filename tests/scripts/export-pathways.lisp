(defun write-to-file (file content)
  "Write a string CONTENT into a file with filename FILE."
  (with-open-file (stream file
                          :direction :output ;; write to disk
                          :if-exists :supersede ;; overwrite
                          :if-does-not-exist :create)
    (write-sequence content stream)))
