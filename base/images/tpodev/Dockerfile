FROM tpobase:temp

# ARGUMENTS
ARG user
ARG userid
ARG passwd

## system-level installs
RUN Rscript -e "install.packages('languageserver')"

## user pass-through
RUN set -e \
    && useradd --no-log-init -m -d /home/$user -u $userid -g users -G sudo $user \
    && echo $user:$passwd | chpasswd

## setup
ADD --chown=$user:users docker_bashrc /home/$user/.bashrc
ADD --chown=$user:users docker_emacs /home/$user/.emacs
RUN su $user -c 'mkdir -p /home/$user/.backup/emacs'
RUN su $user -c 'mkdir -p /home/$user/.local/bin'
RUN su $user -c 'pip3 install moke'
