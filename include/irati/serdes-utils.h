/*
 * Utilities to serialize and deserialize
 *
 *    Eduard Grasa <eduard.grasa@i2cat.net>
 *    Vincenzo Maffione <v.maffione@nextworks.it>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef IRATI_SERDES_H
#define IRATI_SERDES_H

#include "irati/kucommon.h"

#ifdef __cplusplus
extern "C" {
#endif

struct irati_msg_layout {
    unsigned int copylen;
    unsigned int names;
    unsigned int strings;
    unsigned int flow_specs;
    unsigned int dif_configs;
    unsigned int dtp_configs;
    unsigned int dtcp_configs;
    unsigned int buffers;
};

void serialize_string(void **pptr, const char *s);
int deserialize_string(const void **pptr, char **s);

int rina_sername_valid(const char *str);
unsigned rina_name_serlen(const struct rina_name *name);
void serialize_rina_name(void **pptr, const struct name *name);
int deserialize_rina_name(const void **pptr, struct name *name);
void rina_name_free(struct name *name);
void rina_name_move(struct name *dst, struct name *src);
int rina_name_copy(struct name *dst, const struct name *src);
char *rina_name_to_string(const struct name *name);
int rina_name_from_string(const char *str, struct name *name);
int rina_name_cmp(const struct name *one, const struct name *two);
int rina_name_fill(struct name *name, const char *apn,
                   const char *api, const char *aen, const char *aei);
int rina_name_valid(const struct name *name);

int flow_spec_serlen(const struct flow_spec * fspec);
void serialize_flow_spec(void **pptr, const struct flow_spec *fspec);
int deserialize_flow_spec(const void **pptr, struct flow_spec *fspec);
void flow_spec_free(struct flow_spec * fspec);

int policy_parm_serlen(const struct policy_parm * prm);
void serialize_policy_parm(void **pptr, const struct policy_parm *prm);
int deserialize_policy_parm(const void **pptr, struct policy_parm *prm);
void policy_parm_free(struct policy_parm * prm);

int policy_serlen(const struct policy * policy);
void serialize_policy(void **pptr, const struct policy *prm);
int deserialize_policy(const void **pptr, struct policy *prm);
void policy_free(struct policy * policy);

int dtp_config_serlen(const struct dtp_config * dtp_config);
void serialize_dtp_config(void **pptr, const struct dtp_config *dtp_config);
int deserialize_dtp_config(const void **pptr, struct dtp_config *dtp_config);
void dtp_config_free(struct dtp_config * dtp_config);

int window_fctrl_config_serlen(const struct window_fctrl_config * wfc);
void serialize_window_fctrl_config(void **pptr, const struct window_fctrl_config *wfc);
int deserialize_window_fctrl_config(const void **pptr, struct window_fctrl_config *wfc);
void window_fctrl_config_free(struct window_fctrl_config * wfc);

int rate_fctrl_config_serlen(const struct rate_fctrl_config * rfc);
void serialize_rate_fctrl_config(void **pptr, const struct rate_fctrl_config *rfc);
int deserialize_rate_fctrl_config(const void **pptr, struct rate_fctrl_config *rfc);
void rate_fctrl_config_free(struct rate_fctrl_config * rfc);

int dtcp_fctrl_config_serlen(const struct dtcp_fctrl_config * dfc);
void serialize_dtcp_fctrl_config(void **pptr, const struct dtcp_fctrl_config *dfc);
int deserialize_dtcp_fctrl_config(const void **pptr, struct dtcp_fctrl_config *dfc);
void dtcp_fctrl_config_free(struct dtcp_fctrl_config * dfc);

int dtcp_rxctrl_config_serlen(const struct dtcp_rxctrl_config * rxfc);
void serialize_dtcp_rxctrl_config(void **pptr, const struct dtcp_rxctrl_config *rxfc);
int deserialize_dtcp_rxctrl_config(const void **pptr, struct dtcp_rxctrl_config *rxfc);
void dtcp_rxctrl_config_free(struct dtcp_rxctrl_config * rxfc);

int dtcp_config_serlen(const struct dtcp_config * dtcp_config);
void serialize_dctp_config(void **pptr, const struct dtcp_config *dtcp_config);
int deserialize_dtcp_config(const void **pptr, struct dtcp_config *dtcp_config);
void dtcp_config_free(struct dtcp_config * dtcp_config);

int dif_config_serlen(const struct dif_config * dif_config);
void serialize_dif_config(void **pptr, const struct dif_config *dif_config);
int deserialize_dif_config(const void **pptr, struct dif_config *dif_config);

unsigned int irati_msg_serlen(struct irati_msg_layout *numtables,
                              size_t num_entries,
                              const struct irati_msg_base *msg);
unsigned int irati_numtables_max_size(struct irati_msg_layout *numtables,
                                      unsigned int n);
unsigned int serialize_irati_msg(struct irati_msg_layout *numtables,
                                 size_t num_entries,
                                 void *serbuf,
                                 const struct irati_msg_base *msg);
int deserialize_irati_msg(struct irati_msg_layout *numtables, size_t num_entries,
                          const void *serbuf, unsigned int serbuf_len,
                          void *msgbuf, unsigned int msgbuf_len);
void irati_msg_free(struct irati_msg_layout *numtables, size_t num_entries,
                    struct irati_msg_base *msg);

#ifdef __KERNEL__
/* GFP variations of some of the functions above. */
int __rina_name_fill(struct name *name, const char *apn,
                      const char *api, const char *aen, const char *aei,
                      int maysleep);

char * __rina_name_to_string(const struct name *name, int maysleep);
int __rina_name_from_string(const char *str, struct name *name,
                            int maysleep);
#else  /* !__KERNEL__ */
#endif /* !__KERNEL__ */

#ifdef __cplusplus
}
#endif

#endif /* IRATI_SERDES_H */
