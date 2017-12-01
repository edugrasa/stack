/*
 * Reuse Linux Queuing disciplines as RMT scheduling policies
 *
 *    Eduard Grasa <eduard.grasa@i2cat.net>
 *    Miquel Tarzan <miquel.tarzan@i2cat.net>
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

#include <linux/export.h>
#include <linux/module.h>
#include <net/sch_generic.h>
#include <net/pkt_sched.h>
#include <linux/string.h>
#include <linux/version.h>

#define RINA_PREFIX "linux-qdiscs-plugin"

#include "logs.h"
#include "rds/rmem.h"
#include "rmt-ps.h"
#include "policies.h"
#include "debug.h"

#define RINA_LINUX_QDISCS_PS_NAME "linux-qdiscs-ps"
#define FAKE_DEVICE_NAME "fakedev"

struct linux_qdiscs_data {
	struct Qdisc_ops * qdisc_ops;
        unsigned int q_max;
};

struct linux_qdiscs_rmt_queue {
	struct Qdisc * qdisc;
	struct net_device * dev;
};

static struct linux_qdiscs_data * linux_qdiscs_data_create(void)
{
	struct linux_qdiscs_data * tmp;

	tmp = rkzalloc(sizeof(*tmp), GFP_ATOMIC);
	if (!tmp) {
		LOG_ERR("Problems creating Linux Qdiscs data");
		return NULL;
	}

	return tmp;
}

static void linux_qdiscs_data_destroy(struct linux_qdiscs_data * data)
{
	rkfree(data);
}

struct du * linux_qdiscs_dequeue_policy(struct rmt_ps	  *ps,
				   	struct rmt_n1_port *n1_port)
{
	struct linux_qdiscs_rmt_queue *q;
	struct sk_buff * ret;
	struct sock_skb_du * skb_du;

	if (!ps || !n1_port) {
		LOG_ERR("Wrong input parameters");
		return NULL;
	}

	q = n1_port->rmt_ps_queues;
	if (!q) {
		LOG_ERR("Could not find queue for n1_port %u",
			n1_port->port_id);
		return NULL;
	}

	ret = q->qdisc->dequeue(q->qdisc);
	if (!ret) {
		LOG_ERR("Could not dequeue scheduled pdu");
		return NULL;
	}

	skb_du = (struct sock_skb_du *) ret->cb;
	if (!skb_du || !skb_du->du) {
		LOG_ERR("Could not retrieve DU from SKB");
		kfree_skb(ret);
		return NULL;
	}

	return skb_du->du;
}

int linux_qdiscs_enqueue_policy(struct rmt_ps	  *ps,
			   	struct rmt_n1_port *n1_port,
				struct du	  *du)
{
	struct linux_qdiscs_rmt_queue *q;
	int result;
#if LINUX_VERSION_CODE >= KERNEL_VERSION(4,8,0)
	struct sk_buff * to_free = NULL;
#else
#endif

	if (!ps || !n1_port || !du) {
		LOG_ERR("Wrong input parameters");
		return RMT_PS_ENQ_ERR;
	}

	q = n1_port->rmt_ps_queues;
	if (!q) {
		LOG_ERR("Could not find queue for n1_port %u",
			n1_port->port_id);
		du_destroy(du);
		return RMT_PS_ENQ_ERR;
	}

#if LINUX_VERSION_CODE >= KERNEL_VERSION(4,8,0)
	result = q->qdisc->enqueue(du->skb, q->qdisc, &to_free);
#else
	result = q->qdisc->enqueue(du->skb, q->qdisc);
#endif
	if (result ==  NET_XMIT_SUCCESS)
		return RMT_PS_ENQ_SCHED;

	/* Skb has been dropped */
#if LINUX_VERSION_CODE >= KERNEL_VERSION(4,8,0)
#else
	/* Skb has already been freed */
	du->skb = NULL;
#endif

	du_destroy(du);

	return RMT_PS_ENQ_ERR;
}

static void fake_dev_setup(struct net_device *dev)
{
	return;
}

static void * linux_qdiscs_q_create_policy(struct rmt_ps      *ps,
				           struct rmt_n1_port *port)
{
        struct linux_qdiscs_rmt_queue * q;
        struct linux_qdiscs_data * data;

        if (!ps || !port || !ps->priv) {
                LOG_ERR("Wrong input parms for linux_qdiscs_q_create_policy");
		return NULL;
        }

        data = ps->priv;

        q = rkzalloc(sizeof(*q), GFP_ATOMIC);
        if (!q) {
        	LOG_ERR("Problems creating Linux Qdiscs data");
        	return NULL;
        }

        q->dev = alloc_netdev(0, FAKE_DEVICE_NAME,
        		      NET_NAME_UNKNOWN, fake_dev_setup);
        if (!q->dev) {
        	LOG_ERR("Error allocating net device");
        	rkfree(q);
        	return NULL;
        }

	q->qdisc = qdisc_create_dflt(&q->dev->_tx[0], data->qdisc_ops, 0);
	if (!q->qdisc) {
		LOG_ERR("Problems creating qdisc");
		free_netdev(q->dev);
		rkfree(q);
		return NULL;
	}
	q->qdisc->limit = data->q_max;

        return q;
}

static int linux_qdiscs_q_destroy_policy(struct rmt_ps      *ps,
				         struct rmt_n1_port *port)
{
        struct linux_qdiscs_rmt_queue *   q;

        if (!ps || !port || !ps->priv) {
                LOG_ERR("Wrong input parms for rmt_q_destroy_policy");
                return -1;
        }

        q = port->rmt_ps_queues;

        qdisc_destroy(q->qdisc);
	free_netdev(q->dev);
	rkfree(q);

        return 0;
}

static int rmt_config_apply(struct policy_parm * param, void * data)
{
	struct linux_qdiscs_data * tmp;

	tmp = (struct linux_qdiscs_data *) data;

	/* TODO parse RMT policy params and select proper values
	 * for linux_qdiscs_data
	 */

	return 0;
}

static int linux_qdiscs_ps_set_policy_set_param(struct ps_base * bps,
				       	        const char *     name,
						const char *     value)
{
	LOG_ERR("Operation not yet supported");

	return -1;
}

static struct ps_base *
rmt_ps_linux_qdiscs_create(struct rina_component * component)
{
	struct rmt * rmt = rmt_from_component(component);
	struct rmt_ps * ps = rkzalloc(sizeof(*ps), GFP_ATOMIC);
	struct linux_qdiscs_data * data;
	struct rmt_config * rmt_cfg;

	if (!ps)
		return NULL;

	data = linux_qdiscs_data_create();
	if (!data) {
		rkfree(ps);
		return NULL;
	}

	/* FIXME: to be configured using rmt_config */
	data->q_max = 1500;
	data->qdisc_ops = &pfifo_qdisc_ops;

	ps->base.set_policy_set_param = linux_qdiscs_ps_set_policy_set_param;
	ps->dm = rmt;
	ps->priv = data;

	rmt_cfg = rmt_config_get(rmt);
	if (rmt_cfg) {
		policy_for_each(rmt_cfg->policy_set, data, rmt_config_apply);
	} else {
		/* TODO provide a suitable default for all the parameters. */
		LOG_WARN("Missing defaults");
	}

	ps->rmt_dequeue_policy = linux_qdiscs_dequeue_policy;
	ps->rmt_enqueue_policy = linux_qdiscs_enqueue_policy;
	ps->rmt_q_create_policy = linux_qdiscs_q_create_policy;
	ps->rmt_q_destroy_policy = linux_qdiscs_q_destroy_policy;

	LOG_INFO("Loaded Linux Qdiscs policy set and its configuration");

	return &ps->base;
}

static void rmt_ps_linux_qdiscs_destroy(struct ps_base * bps)
{
	struct rmt_ps *ps = container_of(bps, struct rmt_ps, base);

	if (bps) {
		if (ps && ps->priv)
			linux_qdiscs_data_destroy(ps->priv);
	}
}

static struct ps_factory linux_qdiscs_factory = {
		.owner   = THIS_MODULE,
		.create  = rmt_ps_linux_qdiscs_create,
		.destroy = rmt_ps_linux_qdiscs_destroy,
};

static int __init mod_init(void)
{
	int ret;

	strcpy(linux_qdiscs_factory.name, RINA_LINUX_QDISCS_PS_NAME);

	ret = rmt_ps_publish(&linux_qdiscs_factory);
	if (ret) {
		LOG_ERR("Failed to publish policy set factory");
		return -1;
	}

	LOG_INFO("RMT Linux Qdiscs policy set loaded successfully");

	return 0;
}

static void __exit mod_exit(void)
{
	int ret = rmt_ps_unpublish(RINA_LINUX_QDISCS_PS_NAME);

	if (ret) {
		LOG_ERR("Failed to unpublish policy set factory");
		return;
	}

	LOG_INFO("RMT Linux Qdiscs policy set unloaded successfully");
}

module_init(mod_init);
module_exit(mod_exit);

MODULE_DESCRIPTION("RMT Linux Qdiscs policy set");

MODULE_LICENSE("GPL");

MODULE_AUTHOR("Eduard Grasa <eduard.grasa@i2cat.net>");
